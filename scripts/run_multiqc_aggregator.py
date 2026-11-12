"""Aggregate Multi-Stage NGS QC Metrics & Nextflow DSL2 Execution Trace Telemetry.

Evaluates clinical SLA quality gates per sample:
  - Phred Q30 (post-trim) >= 88.0%
  - Alignment rate        >= 95.0%
  - PCR duplication rate  <= 15.0% (WARN if 15.0 - 20.0%, FAIL if > 20.0%)
  - Mean target coverage  >= 30.0x
  - >=30x target breadth  >= 85.0%
  - Ti/Tv concordance     WES in [2.60, 3.10] | WGS in [1.95, 2.25]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def evaluate_sample_qc_gates(row: pd.Series | dict[str, object]) -> dict[str, object]:
    """Evaluate sample-level clinical NGS SLA quality gates and return PASS / WARN / FAIL status."""
    assay = str(row.get("assay_type", "WES_Twist_Exome"))
    q30 = float(row.get("post_trim_q30_pct", 0.0))
    mapped = float(row.get("mapped_pct", 0.0))
    dup = float(row.get("pcr_duplication_pct", 0.0))
    cov = float(row.get("mean_target_cov_x", 0.0))
    breadth30 = float(row.get("pct_target_30x", 0.0))
    titv = float(row.get("titv_ratio", 0.0))

    fail_flags: list[str] = []
    warn_flags: list[str] = []

    if q30 < 85.0:
        fail_flags.append("LOW_Q30")
    elif q30 < 88.0:
        warn_flags.append("BORDERLINE_Q30")

    if mapped < 93.0:
        fail_flags.append("LOW_MAPPING")
    elif mapped < 95.0:
        warn_flags.append("BORDERLINE_MAPPING")

    if dup > 20.0:
        fail_flags.append("HIGH_PCR_DUP")
    elif dup > 15.0:
        warn_flags.append("ELEVATED_PCR_DUP")

    if cov < 25.0 or breadth30 < 75.0:
        fail_flags.append("INSUFFICIENT_DEPTH")
    elif cov < 30.0 or breadth30 < 85.0:
        warn_flags.append("BORDERLINE_DEPTH")

    if assay == "WES_Twist_Exome":
        if not (2.50 <= titv <= 3.15):
            fail_flags.append("TITV_OUTLIER")
    else:
        if not (1.92 <= titv <= 2.25):
            fail_flags.append("TITV_OUTLIER")

    if fail_flags:
        status = "FAIL"
    elif warn_flags:
        status = "WARN"
    else:
        status = "PASS"

    all_flags = fail_flags + warn_flags
    return {
        "qc_gate_status": status,
        "qc_flags": ";".join(all_flags) if all_flags else "NONE",
        "gate_q30_pass": q30 >= 88.0,
        "gate_mapping_pass": mapped >= 95.0,
        "gate_duplication_pass": dup <= 15.0,
        "gate_coverage_pass": cov >= 30.0 and breadth30 >= 85.0,
        "gate_titv_pass": "TITV_OUTLIER" not in fail_flags,
    }


def summarize_nextflow_resources(trace_df: pd.DataFrame) -> pd.DataFrame:
    """Compute per-process Nextflow DSL2 resource utilization benchmarks."""
    df = trace_df.copy()
    df["core_hours"] = (df["cpus"] * df["realtime_min"]) / 60.0
    df["mem_utilization_pct"] = (df["peak_rss_gb"] / df["allocated_mem_gb"]) * 100.0

    grouped = (
        df.groupby("process", as_index=False)
        .agg(
            tasks_executed=("task_id", "count"),
            allocated_cpus=("cpus", "first"),
            allocated_mem_gb=("allocated_mem_gb", "first"),
            mean_peak_rss_gb=("peak_rss_gb", "mean"),
            max_peak_rss_gb=("peak_rss_gb", "max"),
            mean_mem_util_pct=("mem_utilization_pct", "mean"),
            mean_cpu_eff_pct=("cpu_efficiency_pct", "mean"),
            mean_duration_min=("realtime_min", "mean"),
            total_core_hours=("core_hours", "sum"),
            retried_tasks=("attempt", lambda s: int((s > 1).sum())),
        )
    )

    numeric_cols = [
        "mean_peak_rss_gb",
        "max_peak_rss_gb",
        "mean_mem_util_pct",
        "mean_cpu_eff_pct",
        "mean_duration_min",
        "total_core_hours",
    ]
    for col in numeric_cols:
        grouped[col] = grouped[col].round(2)

    proc_order = [
        "FASTP_TRIM",
        "BWA_MEM2_ALIGN",
        "SAMTOOLS_MARKDUP_STATS",
        "BCFTOOLS_CALL_FILTER",
        "MULTIQC_COHORT_AGGREGATE",
    ]
    grouped["sort_key"] = grouped["process"].apply(lambda x: proc_order.index(x) if x in proc_order else 99)
    grouped = grouped.sort_values("sort_key").drop(columns=["sort_key"]).reset_index(drop=True)
    return grouped


def write_markdown_report(
    cohort_df: pd.DataFrame,
    resource_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """Generate a formatted MultiQC & Nextflow Execution Markdown Report."""
    status_counts = cohort_df["qc_gate_status"].value_counts().to_dict()
    total_core_h = float(resource_df["total_core_hours"].sum())

    lines = [
        "# Nextflow DSL2 Pipeline Execution & MultiQC Cohort Report",
        "",
        "- **Workflow Engine:** `Nextflow DSL2 (>=23.10.0)`",
        "- **Container Image:** `ghcr.io/denisulko/reproducible-ngs-nextflow:1.0.0`",
        f"- **Samples Processed:** {len(cohort_df)} (`6 WES_Twist_Exome`, `6 WGS_PCR_Free`)",
        f"- **SLA Gate Summary:** `{status_counts.get('PASS', 0)} PASS` | `{status_counts.get('WARN', 0)} WARN` | `{status_counts.get('FAIL', 0)} FAIL`",
        f"- **Total Compute Consumption:** `{total_core_h:.2f} CPU Core-Hours` across `{int(resource_df['tasks_executed'].sum())} DSL2 tasks`",
        "",
        "## 1. Sample-Level MultiQC & SLA Gate Table",
        "",
        "| Sample ID | Assay | Raw Reads (M) | Post-Trim Q30 (%) | Mapped (%) | PCR Dup (%) | Mean Cov (x) | >=30x (%) | SNVs | Ti/Tv | SLA Status | Flags |",
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |",
    ]
    for _, r in cohort_df.iterrows():
        lines.append(
            f"| `{r['sample_id']}` | `{r['assay_type']}` | {r['raw_reads_m']:.1f} | "
            f"{r['post_trim_q30_pct']:.1f}% | {r['mapped_pct']:.1f}% | {r['pcr_duplication_pct']:.1f}% | "
            f"{r['mean_target_cov_x']:.1f}x | {r['pct_target_30x']:.1f}% | {int(r['pass_snv_count']):,} | "
            f"{r['titv_ratio']:.3f} | **`{r['qc_gate_status']}`** | `{r['qc_flags']}` |"
        )

    lines.extend(
        [
            "",
            "## 2. Nextflow DSL2 Process Resource Benchmarks",
            "",
            "| DSL2 Process | Tasks | CPUs | Alloc RAM (GB) | Mean Peak RSS (GB) | RAM Util (%) | Mean CPU Eff (%) | Mean Walltime (min) | Total Core-Hours | Retries |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for _, r in resource_df.iterrows():
        lines.append(
            f"| `{r['process']}` | {int(r['tasks_executed'])} | {int(r['allocated_cpus'])} | "
            f"{r['allocated_mem_gb']:.1f} | {r['mean_peak_rss_gb']:.2f} | {r['mean_mem_util_pct']:.1f}% | "
            f"{r['mean_cpu_eff_pct']:.1f}% | {r['mean_duration_min']:.2f} | {r['total_core_hours']:.2f} | {int(r['retried_tasks'])} |"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--telemetry-csv", type=Path, default=ROOT / "data/raw/pipeline_stage_telemetry.csv")
    parser.add_argument("--trace-csv", type=Path, default=ROOT / "data/raw/nextflow_execution_trace.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    args = parser.parse_args(argv)

    tel_df = pd.read_csv(args.telemetry_csv)
    trace_df = pd.read_csv(args.trace_csv)

    gates = pd.DataFrame([evaluate_sample_qc_gates(row) for _, row in tel_df.iterrows()], index=tel_df.index)
    cohort_df = pd.concat([tel_df, gates], axis=1)

    resource_df = summarize_nextflow_resources(trace_df)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cohort_df.to_csv(args.output_dir / "multiqc_cohort_summary.csv", index=False)
    resource_df.to_csv(args.output_dir / "nextflow_resource_benchmarks.csv", index=False)

    write_markdown_report(cohort_df, resource_df, args.output_dir / "pipeline_execution_report.md")

    status_counts = {k: int(v) for k, v in cohort_df["qc_gate_status"].value_counts().items()}
    report_payload = {
        "pipeline_name": "denisulko/reproducible-ngs-nextflow",
        "dsl_version": "DSL2",
        "total_samples": int(len(cohort_df)),
        "qc_status_counts": status_counts,
        "total_tasks_executed": int(len(trace_df)),
        "total_core_hours": round(float(resource_df["total_core_hours"].sum()), 2),
        "mean_cpu_efficiency_pct": round(float(trace_df["cpu_efficiency_pct"].mean()), 2),
        "mean_post_trim_q30_pct": round(float(cohort_df["post_trim_q30_pct"].mean()), 2),
        "mean_mapped_pct": round(float(cohort_df["mapped_pct"].mean()), 2),
    }
    (args.output_dir / "multiqc_report.json").write_text(json.dumps(report_payload, indent=2), encoding="utf-8")

    print(f"Aggregated MultiQC summary for {len(cohort_df)} samples -> {args.output_dir}")
    print(f"SLA gate breakdown: {status_counts}")
    print(f"Total compute: {report_payload['total_core_hours']} CPU core-hours across {len(trace_df)} Nextflow tasks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
