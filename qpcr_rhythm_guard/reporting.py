"""Inspectable CSV/JSON/HTML and Matplotlib figures, without external resources."""
import csv
import hashlib
import html
import json
import math
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

REASONS = {
    "design_unknown_descriptive_only": "重复设计未确认，仅作描述",
    "design_repeated_descriptive_only": "重复测量设计，仅作描述",
    "fewer_than_four_distinct_phases": "有效相位少于 4 个",
    "time_span_below_three_quarters_period": "有效时间跨度不足 0.75 个周期",
    "fewer_than_two_subjects_at_a_time": "至少一个输入时间点的有效独立样本不足 2 个",
    "normalization_quality_gate_failed": "技术重复或归一化质控未通过推断门槛",
    "perfect_fit_no_residual_variance": "残差方差近乎为零，未计算 P 值",
    "multiple_runs_without_validated_batch_model": "存在多个批次，首版未进行批次建模",
    "subject_repeated_across_time": "同一生物学单位出现在多个时间点",
    "no_complete_observations": "没有完整归一化观测",
    "rank_deficient_design": "采样设计不能唯一估计模型参数",
    "ill_conditioned_design": "采样设计数值不稳定",
    "insufficient_residual_degrees_of_freedom": "残差自由度不足",
}

QC_TEXT = {
    "missing_ct": "输入包含缺失或未检出 Ct，原值已保留，未补为零。",
    "normalization_assumptions": "单内参稳定性和扩增效率假设需另行验证。",
    "missing_observations_excluded": "缺少完整归一化值的记录没有参与拟合；缺失可能并非随机。",
    "inferential_test_withheld": "当前模型仅作描述性拟合，未进行推断性检验。",
    "sample_level_analysis": "技术重复先在样本内聚合；有效 N 计样本观测，独立单位数另行保留。",
    "high_technical_sd": "技术重复 Ct 的标准差超过设定阈值，保留描述值并暂停相关推断。",
    "single_valid_technical_replicate": "仅剩一个有效技术重复，无法估计技术标准差。",
    "reference_gene_absent": "未找到指定内参，无法计算归一化值。",
    "no_target_genes": "排除内参后没有目标基因。",
    "missing_calibrator": "指定校准组和时间没有完整的目标基因与内参配对。",
    "incomplete_calibrators_excluded": "不完整校准样本未计入校准均值，缺失记录仍保留。",
    "multiple_runs": "输入包含多个批次；首版没有批次校正，所有模型仅作描述。",
    "repeated_subjects": "存在跨时间重复取样，暂停普通独立样本 F 检验。",
    "rq_out_of_range": "RQ 超出可表示范围；保留 log2RQ 描述值并暂停相关推断。",
}


def scalar(value):
    if isinstance(value, np.generic):
        return scalar(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): scalar(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [scalar(v) for v in value]
    return value


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row)) or ["message"]
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = {}
            for key, value in row.items():
                if isinstance(value, (dict, list, tuple)):
                    value = json.dumps(value, ensure_ascii=False)
                # Labels are never executed as spreadsheet formulas. JSON retains exact strings.
                if isinstance(value, str) and value.startswith(("=", "+", "-", "@", "\t", "\r")):
                    value = "'" + value
                out[key] = value
            writer.writerow(out)


def display(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.5g}"
    if isinstance(value, (dict, list)):
        value = json.dumps(value, ensure_ascii=False)
    return html.escape(str(value))


def table(rows, fields):
    header = "".join(f"<th>{html.escape(label)}</th>" for key, label in fields)
    body = "".join("<tr>" + "".join(f"<td>{display(row.get(key))}</td>" for key, _ in fields) + "</tr>" for row in rows)
    return f'<div class="table-wrap"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>'


def figures(data, out):
    files = []
    colors = ["#087e8b", "#cd5c35", "#7562a2", "#486a3e"]
    for index, gene in enumerate(sorted({r["gene"] for r in data["fits"]}), 1):
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), layout="constrained")
        subset = [r for r in data["normalized"] if r["gene"] == gene and r.get("log2rq") is not None]
        for j, group in enumerate(sorted({r["group"] for r in subset})):
            rows = [r for r in subset if r["group"] == group]
            color = colors[j % len(colors)]
            x = np.array([r["time_h"] for r in rows])
            y = np.array([r["log2rq"] for r in rows])
            axes[0].scatter(x, [r["rq"] for r in rows], color=color, alpha=.7, s=19, label=group)
            axes[1].scatter(x, y, color=color, alpha=.7, s=19, label=group)
            fit = next(r for r in data["fits"] if r["gene"] == gene and r["group"] == group)
            if fit.get("mesor") is not None:
                xx = np.linspace(x.min(), x.max(), 240)
                yy = fit["mesor"] + fit["cos_coef"] * np.cos(2 * np.pi * xx / fit["period_h"]) + fit["sin_coef"] * np.sin(2 * np.pi * xx / fit["period_h"])
                axes[1].plot(xx, yy, color=color, linewidth=2)
        for ax, title, ylabel in zip(axes, ["Observed normalized values", "Fixed-period descriptive fit"], ["RQ", "log2(RQ)"]):
            ax.set(title=title, xlabel="Elapsed time (h)", ylabel=ylabel)
            ax.grid(alpha=.16)
            ax.spines[["top", "right"]].set_visible(False)
            if subset:
                ax.legend(frameon=False, fontsize=8)
        fig.suptitle(f"{gene} | points = biological samples; no extrapolation", fontsize=11)
        name = f"figure_{index:03d}"
        fig.savefig(out / f"{name}.svg")
        fig.savefig(out / f"{name}.png", dpi=160)
        plt.close(fig)
        files.append((gene, f"{name}.svg"))
    return files


def export_results(result, out, source, config_path):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    data = scalar(result)
    for name in ("normalized", "fits", "qc"):
        write_csv(out / f"{name}.csv", data[name])
    (out / "results.json").write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    plots = figures(data, out)
    valid = sum(r.get("log2rq") is not None for r in data["normalized"])
    tested = sum(r.get("p_value") is not None for r in data["fits"])
    fit_rows = []
    for row in data["fits"]:
        row = dict(row)
        row["status"] = {"tested": "已检验", "descriptive": "仅描述", "not_testable": "不可检验"}.get(row["status"], row["status"])
        if row.get("reason"):
            row["reason"] = "；".join(REASONS.get(reason, reason) for reason in row["reason"].split("; "))
        fit_rows.append(row)
    fit_table = table(fit_rows, [("gene", "基因"), ("group", "组别"), ("status", "状态"), ("n_input", "输入 N"), ("n", "有效 N"), ("amplitude", "振幅 / log2RQ"), ("peak_h", "峰相位 / h"), ("p_value", "原始 P"), ("q_value", "BH q"), ("reason", "限制或说明")])
    qc_rows = []
    for row in data["qc"]:
        row = dict(row)
        row["context"] = " · ".join(f"{key}={value}" for key, value in row.items() if key not in ("level", "code", "message"))
        row["message"] = QC_TEXT.get(row["code"], row["message"])
        row["level"] = {"warning": "需核对", "info": "说明", "error": "错误"}.get(row["level"], row["level"])
        qc_rows.append(row)
    qc_table = table(qc_rows, [("level", "级别"), ("message", "说明"), ("context", "对应记录")])
    charts = "".join(f'<section><h2>{html.escape(gene)}</h2><img src="{file}" alt="{html.escape(gene)} 观测值与固定周期拟合图"></section>' for gene, file in plots)
    config_html = html.escape(json.dumps(data["config"], ensure_ascii=False, indent=2))
    report = f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>qPCR Rhythm Guard · 分析报告</title>
<style>body{{margin:0;background:#f2f6f5;color:#203733;font:15px/1.65 "Segoe UI","Microsoft YaHei",sans-serif}}main{{max-width:1180px;margin:auto;padding:42px 24px}}header{{border-top:5px solid #087e8b;padding:24px 0}}.eyebrow{{letter-spacing:2px;color:#087e8b;font-size:12px}}h1{{font-size:36px;line-height:1.2;margin:12px 0}}h2{{font-size:21px;margin-top:0}}.muted{{color:#627670}}.cards{{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:24px 0}}.card,section{{background:#fff;border:1px solid #d9e4df;border-radius:10px;padding:24px;margin-bottom:20px}}.card strong{{display:block;font-size:30px;color:#087e8b}}.note{{border-left:4px solid #cb8147;background:#fff5e9;padding:16px 20px}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{text-align:left;padding:10px 12px;border-bottom:1px solid #e2e9e5;vertical-align:top}}th{{white-space:nowrap;background:#edf4f1}}td{{overflow-wrap:anywhere}}.table-wrap{{overflow-x:auto}}img{{width:100%;height:auto}}a{{color:#086d79}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}}footer{{font-size:12px;color:#627670}}@media(max-width:700px){{.cards{{grid-template-columns:1fr}}main{{padding:20px 12px}}h1{{font-size:28px}}section{{padding:16px}}}}</style>
<main><header><div class="eyebrow">LOCAL ANALYSIS / v0.1.0</div><h1>qPCR Rhythm Guard</h1><p class="muted">从原始 Ct 到可追溯的质控与节律拟合。数据在本机处理。</p><p>输入：{html.escape(Path(source).name)} · 内参：{html.escape(str(data['config']['reference_gene']))} · 固定周期：{display(data['config']['period_h'])} h</p></header>
<div class="note">这是研究辅助结果。固定周期拟合不能证明内源性昼夜节律。NA 表示未计算或不可检验，不能解释为无变化。经典 F 检验的残差假设仍需人工评估。</div>
<div class="cards"><div class="card"><strong>{valid}</strong>有效归一化记录（跨基因）</div><div class="card"><strong>{tested} / {len(data['fits'])}</strong>可检验 / 全部模型</div><div class="card"><strong>{len(data['qc'])}</strong>质控提示</div></div>
<section><h2>模型与检验</h2><p>每个点代表一个生物学样本；技术重复已聚合。振幅在 log2RQ 尺度，峰相位相对于输入时间 0。BH 家族为本次运行全部可检验的基因 × 组别，不能用两组 P 值的显著与否判断组间差异。</p>{fit_table}</section>
{charts}<section><h2>质控记录</h2><p>缺失和未检出没有补成零。提示需要结合原始扩增曲线、内参稳定性和实验记录核对。</p>{qc_table}</section>
<section><h2>下载与复核</h2><p><a href="normalized.csv">归一化数据 CSV</a> · <a href="fits.csv">模型结果 CSV</a> · <a href="qc.csv">质控记录 CSV</a> · <a href="results.json">完整数据 JSON</a> · <a href="manifest.json">参数与文件校验值</a></p><details><summary>运行参数</summary><pre>{config_html}</pre></details></section>
<footer>qPCR Rhythm Guard 0.1.0 · 研究原型 · 等效率转换与单内参稳定性需要实验验证；本软件不认证 MIQE 合规。示例文件 synthetic.csv 全部为人工生成数据。</footer></main></html>'''
    (out / "report.html").write_text(report, encoding="utf-8")
    manifest = {
        "project": "qPCR Rhythm Guard", "version": "0.1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_name": Path(source).name, "source_sha256": hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(Path(config_path).read_bytes()).hexdigest(),
        "config": data["config"], "python": platform.python_version(),
        "dependencies": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "bh_family": "All finite p-values from all gene-by-group models in this invocation",
        "n_planned_models": len(data["fits"]), "n_tested_models": tested,
        "csv_note": "Text starting with formula characters is prefixed by apostrophe; JSON preserves identifiers.",
        "files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.is_file() and p.name != "manifest.json"},
    }
    (out / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
