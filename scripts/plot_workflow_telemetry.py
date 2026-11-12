"""Generate 300 DPI MultiQC Cohort Telemetry & Nextflow DSL2 Resource Profiling Figures.

Outputs in results/:
  1. multiqc_cohort_dashboard.png   — 4-panel MultiQC sample retention, coverage, Q30 & SLA gate matrix
  2. nextflow_resource_profiling.png — 2-panel DSL2 Peak RSS memory utilization & CPU/walltime telemetry
  3. variant_concordance_titv.png   — 2-panel WES/WGS variant yield & Ti/Tv vs. Het/Hom concordance
"""

from __future__ import annotations

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

STATUS_COLORS = {
    "PASS": "#16a34a",
    "WARN": "#f59e0b",
    "FAIL": "#dc2626",
}

PROCESS_COLORS = {
    "FASTP_TRIM": "#0284c7",
    "BWA_MEM2_ALIGN": "#dc2626",
    "SAMTOOLS_MARKDUP_STATS": "#7c3aed",
    "BCFTOOLS_CALL_FILTER": "#ea580c",
    "MULTIQC_COHORT_AGGREGATE": "#16a34a",
}


def plot_multiqc_dashboard(cohort_df: pd.DataFrame, output_path: Path) -> None:
    """Create 4-panel MultiQC cohort quality control and SLA gate figure."""
    fig, axes = plt.subplots(2, 2, figsize=(16.0, 11.5), dpi=300)
    ax1, ax2, ax3, ax4 = axes.flatten()

    samples = cohort_df["sample_id"].tolist()
    x = np.arange(len(samples))

    # Panel A: Read Retention % Across Pipeline Stages
    trim_pct = (cohort_df["trimmed_reads_m"] / cohort_df["raw_reads_m"]) * 100.0
    map_pct = (cohort_df["mapped_reads_m"] / cohort_df["raw_reads_m"]) * 100.0
    dedup_pct = (cohort_df["dedup_reads_m"] / cohort_df["raw_reads_m"]) * 100.0

    w = 0.26
    ax1.bar(x - w, trim_pct, width=w, color="#38bdf8", edgecolor="#0f172a", label="Post-fastp Trimmed (%)")
    ax1.bar(x, map_pct, width=w, color="#2563eb", edgecolor="#0f172a", label="BWA-MEM2 Mapped (%)")
    ax1.bar(x + w, dedup_pct, width=w, color="#1e3a8a", edgecolor="#0f172a", label="Unique Deduplicated (%)")
    ax1.axhline(80.0, color="#dc2626", linestyle="--", lw=1.1, label="Min Dedup Yield Gate (80%)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(samples, rotation=35, ha="right", fontsize=9.5, fontweight="bold")
    ax1.set_ylim(50, 106)
    ax1.set_ylabel("% of Raw Sequenced Read Pairs", fontsize=10.5, fontweight="bold")
    ax1.set_title("A. Multi-Stage Read Retention Funnel (% of Raw Reads)", fontsize=12, fontweight="bold", pad=10)
    ax1.legend(loc="lower left", fontsize=8.5, frameon=True, facecolor="white")
    ax1.grid(axis="y", linestyle="--", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Mean Target Coverage vs. >=30x Breadth
    for status in ["PASS", "WARN", "FAIL"]:
        sub = cohort_df[cohort_df["qc_gate_status"] == status]
        if sub.empty:
            continue
        for assay, marker in [("WES_Twist_Exome", "o"), ("WGS_PCR_Free", "s")]:
            sub_a = sub[sub["assay_type"] == assay]
            if sub_a.empty:
                continue
            ax2.scatter(
                sub_a["mean_target_cov_x"],
                sub_a["pct_target_30x"],
                color=STATUS_COLORS[status],
                marker=marker,
                s=85,
                alpha=0.92,
                edgecolors="#0f172a",
                linewidths=0.8,
                label=f"{status} - {assay.split('_')[0]}",
            )

    for _, r in cohort_df[cohort_df["qc_gate_status"] != "PASS"].iterrows():
        ax2.annotate(
            f"{r['sample_id']} ({r['qc_gate_status']})",
            xy=(float(r["mean_target_cov_x"]), float(r["pct_target_30x"])),
            xytext=(10, -8 if r["sample_id"] == "NGS_S08" else 8),
            textcoords="offset points",
            fontsize=8.8,
            fontweight="bold",
            color="#991b1b" if r["qc_gate_status"] == "FAIL" else "#92400e",
            bbox=dict(boxstyle="round,pad=0.2", fc="#fef2f2", ec="#f87171", lw=0.7),
            arrowprops=dict(arrowstyle="->", color="#991b1b", lw=0.7),
        )

    ax2.axvline(30.0, color="#dc2626", linestyle="--", lw=1.1, label="Min Mean Depth (30x)")
    ax2.axhline(85.0, color="#ea580c", linestyle=":", lw=1.1, label="Min >=30x Breadth (85%)")
    ax2.set_xlabel("Mean Target Sequencing Depth (x)", fontsize=10.5, fontweight="bold")
    ax2.set_ylabel("Target Bases >= 30x Coverage (%)", fontsize=10.5, fontweight="bold")
    ax2.set_title("B. Target Coverage Depth vs. >=30x Completeness", fontsize=12, fontweight="bold", pad=10)
    ax2.legend(loc="lower right", fontsize=8.2, frameon=True, facecolor="white")
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Panel C: Pre-Trim vs. Post-Trim Phred Q30 Dumbbell Plot
    for idx, r in cohort_df.iterrows():
        ax3.plot(
            [idx, idx],
            [r["pre_trim_q30_pct"], r["post_trim_q30_pct"]],
            color="#64748b",
            lw=1.8,
            zorder=1,
        )
    ax3.scatter(
        x,
        cohort_df["pre_trim_q30_pct"],
        color="#94a3b8",
        edgecolors="#0f172a",
        s=65,
        zorder=2,
        label="Pre-Trim Raw Q30 (%)",
    )
    ax3.scatter(
        x,
        cohort_df["post_trim_q30_pct"],
        color=[STATUS_COLORS[s] for s in cohort_df["qc_gate_status"]],
        edgecolors="#0f172a",
        s=78,
        zorder=3,
        label="Post-fastp Q30 (%) [Colored by SLA]",
    )
    ax3.axhline(88.0, color="#dc2626", linestyle="--", lw=1.2, label="Clinical Q30 Threshold (88%)")
    ax3.set_xticks(x)
    ax3.set_xticklabels(samples, rotation=35, ha="right", fontsize=9.5, fontweight="bold")
    ax3.set_ylim(73, 100)
    ax3.set_ylabel("Phred Quality >= Q30 Bases (%)", fontsize=10.5, fontweight="bold")
    ax3.set_title("C. fastp Base Quality Gain (Pre- vs. Post-Trim Q30)", fontsize=12, fontweight="bold", pad=10)
    ax3.legend(loc="lower left", fontsize=8.5, frameon=True, facecolor="white")
    ax3.grid(True, linestyle="--", alpha=0.35)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)

    # Panel D: Clinical SLA Quality Gate Matrix (12 Samples x 5 Gates)
    gate_cols = [
        "gate_q30_pass",
        "gate_mapping_pass",
        "gate_duplication_pass",
        "gate_coverage_pass",
        "gate_titv_pass",
    ]
    gate_labels = ["Q30 >= 88%", "Mapped >= 95%", "PCR Dup <= 15%", "Depth >= 30x", "Ti/Tv Range"]
    mat = cohort_df[gate_cols].astype(int).values.T  # shape: (5, 12)

    cmap = ListedColormap(["#ef4444", "#22c55e"])
    ax4.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax4.set_xticks(x)
    ax4.set_xticklabels(samples, rotation=35, ha="right", fontsize=9.5, fontweight="bold")
    ax4.set_yticks(np.arange(len(gate_labels)))
    ax4.set_yticklabels(gate_labels, fontsize=10, fontweight="bold")
    ax4.set_title("D. Automated Sample SLA Gate Matrix (PASS / FAIL)", fontsize=12, fontweight="bold", pad=10)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            txt = "PASS" if mat[i, j] == 1 else "FAIL"
            ax4.text(
                j,
                i,
                txt,
                ha="center",
                va="center",
                fontsize=8.2,
                fontweight="bold",
                color="white",
            )

    fig.tight_layout(pad=2.2)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_nextflow_resources(trace_df: pd.DataFrame, resource_df: pd.DataFrame, output_path: Path) -> None:
    """Plot 2-panel Nextflow DSL2 Peak RSS Memory Utilization and Walltime vs. CPU Efficiency."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.6, 6.2), dpi=300)

    # Panel A: Peak RSS vs. Allocated RAM per DSL2 Process
    procs = resource_df["process"].tolist()
    y_pos = np.arange(len(procs))
    alloc = resource_df["allocated_mem_gb"].values
    mean_rss = resource_df["mean_peak_rss_gb"].values
    max_rss = resource_df["max_peak_rss_gb"].values

    ax1.barh(y_pos, alloc, color="#e2e8f0", edgecolor="#64748b", height=0.62, label="Allocated Container Memory Limit (GB)")
    ax1.barh(
        y_pos,
        mean_rss,
        color=[PROCESS_COLORS.get(p, "#2563eb") for p in procs],
        edgecolor="#0f172a",
        height=0.44,
        alpha=0.92,
        label="Mean Peak RSS Used (GB)",
    )
    ax1.scatter(max_rss, y_pos, color="#0f172a", marker="D", s=45, zorder=4, label="Max Peak RSS Observed (GB)")

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(procs, fontsize=9.8, fontweight="bold")
    ax1.invert_yaxis()
    ax1.set_xlabel("Memory (GB)", fontsize=11, fontweight="bold")
    ax1.set_title("A. Nextflow DSL2 Process Memory Footprint vs. Reservation", fontsize=12.2, fontweight="bold", pad=12)
    ax1.set_xlim(0, 28)
    ax1.legend(loc="lower right", fontsize=8.5, frameon=True, facecolor="white")
    ax1.grid(axis="x", linestyle="--", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    for idx, r in resource_df.iterrows():
        ax1.text(
            float(r["allocated_mem_gb"]) + 0.4,
            idx,
            f"{r['mean_mem_util_pct']:.1f}% util",
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold",
            color="#0f172a",
        )

    # Panel B: Task Walltime vs. CPU Efficiency Scatter across 49 Tasks
    for proc in procs:
        sub = trace_df[trace_df["process"] == proc]
        if sub.empty:
            continue
        ax2.scatter(
            sub["realtime_min"],
            sub["cpu_efficiency_pct"],
            color=PROCESS_COLORS.get(proc, "#2563eb"),
            s=58,
            alpha=0.88,
            edgecolors="#0f172a",
            linewidths=0.6,
            label=f"{proc} (n={len(sub)})",
        )

    # Highlight retried task (attempt == 2)
    retried = trace_df[trace_df["attempt"] > 1]
    for _, r in retried.iterrows():
        ax2.annotate(
            f"Auto-Retry (attempt={int(r['attempt'])})\n{r['process']} [{r['tag']}]",
            xy=(float(r["realtime_min"]), float(r["cpu_efficiency_pct"])),
            xytext=(-125, -28),
            textcoords="offset points",
            fontsize=8.3,
            fontweight="bold",
            color="#7f1d1d",
            bbox=dict(boxstyle="round,pad=0.22", fc="#fef2f2", ec="#f87171", lw=0.7),
            arrowprops=dict(arrowstyle="->", color="#991b1b", lw=0.8),
        )

    ax2.axhline(80.0, color="#16a34a", linestyle="--", lw=1.1, label="Target CPU Efficiency (>=80%)")
    ax2.set_xlabel("Task Execution Walltime (Minutes)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("CPU Core Utilization Efficiency (%)", fontsize=11, fontweight="bold")
    ax2.set_title("B. Nextflow Task Execution Walltime vs. CPU Efficiency", fontsize=12.2, fontweight="bold", pad=12)
    ax2.set_ylim(75, 102)
    ax2.legend(loc="lower right", fontsize=8.2, frameon=True, facecolor="white")
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout(pad=2.0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_variant_concordance(cohort_df: pd.DataFrame, output_path: Path) -> None:
    """Plot 2-panel Variant Yield by Assay Type and Ti/Tv vs. Het/Hom Concordance."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.6, 6.2), dpi=300)

    # Panel A: Normalized Variant Yield (log10 SNVs and Indels per Sample)
    samples = cohort_df["sample_id"].tolist()
    x = np.arange(len(samples))
    snv_log = np.log10(cohort_df["pass_snv_count"].clip(lower=1))
    indel_log = np.log10(cohort_df["pass_indel_count"].clip(lower=1))

    w = 0.36
    ax1.bar(x - w / 2, snv_log, width=w, color="#2563eb", edgecolor="#0f172a", label="PASS SNVs (log10 count)")
    ax1.bar(x + w / 2, indel_log, width=w, color="#f97316", edgecolor="#0f172a", label="PASS Indels (log10 count)")
    ax1.axvline(5.5, color="#475569", linestyle="--", lw=1.2)
    ax1.text(2.5, 6.72, "WES Twist Exome\n(~41k SNVs / sample)", ha="center", fontsize=9.5, fontweight="bold", color="#1e3a8a")
    ax1.text(8.5, 6.72, "WGS PCR-Free\n(~3.8M SNVs / sample)", ha="center", fontsize=9.5, fontweight="bold", color="#1e3a8a")

    ax1.set_xticks(x)
    ax1.set_xticklabels(samples, rotation=35, ha="right", fontsize=9.5, fontweight="bold")
    ax1.set_ylim(3.2, 7.2)
    ax1.set_ylabel("log10(Variant Calls Passing QUAL>=30, DP>=15)", fontsize=10.5, fontweight="bold")
    ax1.set_title("A. Germline SNV & Indel Yield Across WES and WGS Cohorts", fontsize=12.2, fontweight="bold", pad=12)
    ax1.legend(loc="lower right", fontsize=8.8, frameon=True, facecolor="white")
    ax1.grid(axis="y", linestyle="--", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Ti/Tv Ratio vs. Het/Hom Ratio with Benchmark Windows
    ax2.axvspan(2.00, 2.15, color="#dbeafe", alpha=0.55, label="Expected WGS Ti/Tv Window (2.00 - 2.15)")
    ax2.axvspan(2.70, 3.00, color="#fef3c7", alpha=0.55, label="Expected WES Coding Ti/Tv Window (2.70 - 3.00)")

    for assay, marker in [("WES_Twist_Exome", "o"), ("WGS_PCR_Free", "s")]:
        sub = cohort_df[cohort_df["assay_type"] == assay]
        for _, r in sub.iterrows():
            st = str(r["qc_gate_status"])
            ax2.scatter(
                r["titv_ratio"],
                r["het_hom_ratio"],
                color=STATUS_COLORS[st],
                marker=marker,
                s=90,
                edgecolors="#0f172a",
                linewidths=0.8,
                zorder=3,
            )

    # Proxy legend items for assay & status
    ax2.scatter([], [], color="#16a34a", marker="o", s=75, edgecolors="#0f172a", label="WES Twist Exome (PASS)")
    ax2.scatter([], [], color="#16a34a", marker="s", s=75, edgecolors="#0f172a", label="WGS PCR-Free (PASS)")
    ax2.scatter([], [], color="#f59e0b", marker="s", s=75, edgecolors="#0f172a", label="WARN Sample (NGS_S08)")
    ax2.scatter([], [], color="#dc2626", marker="s", s=75, edgecolors="#0f172a", label="FAIL Outlier (NGS_S11)")

    outlier = cohort_df[cohort_df["sample_id"] == "NGS_S11"].iloc[0]
    ax2.annotate(
        f"NGS_S11 (Ti/Tv={outlier['titv_ratio']:.2f})\nLow-Depth / FP Drift",
        xy=(float(outlier["titv_ratio"]), float(outlier["het_hom_ratio"])),
        xytext=(12, -22),
        textcoords="offset points",
        fontsize=8.5,
        fontweight="bold",
        color="#7f1d1d",
        bbox=dict(boxstyle="round,pad=0.22", fc="#fef2f2", ec="#f87171", lw=0.7),
        arrowprops=dict(arrowstyle="->", color="#991b1b", lw=0.8),
    )

    ax2.set_xlabel("Transition / Transversion Ratio (Ti/Tv)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Heterozygous / Homozygous Alternative Ratio (Het/Hom)", fontsize=10.5, fontweight="bold")
    ax2.set_title("B. Variant Call Set Concordance (Ti/Tv vs. Het/Hom Ratio)", fontsize=12.2, fontweight="bold", pad=12)
    ax2.set_xlim(1.75, 3.12)
    ax2.set_ylim(1.42, 1.85)
    ax2.legend(loc="upper left", fontsize=8.2, frameon=True, facecolor="white")
    ax2.grid(True, linestyle="--", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout(pad=2.0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-csv", type=Path, default=ROOT / "results/multiqc_cohort_summary.csv")
    parser.add_argument("--resource-csv", type=Path, default=ROOT / "results/nextflow_resource_benchmarks.csv")
    parser.add_argument("--trace-csv", type=Path, default=ROOT / "data/raw/nextflow_execution_trace.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    args = parser.parse_args(argv)

    cohort_df = pd.read_csv(args.cohort_csv)
    resource_df = pd.read_csv(args.resource_csv)
    trace_df = pd.read_csv(args.trace_csv)

    p1 = args.output_dir / "multiqc_cohort_dashboard.png"
    p2 = args.output_dir / "nextflow_resource_profiling.png"
    p3 = args.output_dir / "variant_concordance_titv.png"

    plot_multiqc_dashboard(cohort_df, p1)
    plot_nextflow_resources(trace_df, resource_df, p2)
    plot_variant_concordance(cohort_df, p3)

    print(f"Saved 300 DPI plots:\n  - {p1}\n  - {p2}\n  - {p3}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
