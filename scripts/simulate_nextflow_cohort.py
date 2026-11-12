"""Generate a 12-Sample Clinical NGS Cohort Samplesheet, Multi-Stage Telemetry & Nextflow Trace.

Produces:
  - data/raw/samplesheet.csv              (Nextflow DSL2 input channel manifest)
  - data/raw/pipeline_stage_telemetry.csv (fastp, BWA-MEM2, samtools markdup, bcftools metrics)
  - data/raw/nextflow_execution_trace.csv (49 Nextflow task execution records across 5 DSL2 processes)
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def generate_nextflow_cohort(output_dir: Path, seed: int = 42) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Generate deterministic 12-sample NGS cohort, stage telemetry, and Nextflow execution trace."""
    rng = random.Random(seed)
    output_dir.mkdir(parents=True, exist_ok=True)

    samplesheet_rows: list[dict[str, object]] = []
    telemetry_rows: list[dict[str, object]] = []
    trace_rows: list[dict[str, object]] = []

    task_id = 1
    for idx in range(1, 13):
        sid = f"NGS_S{idx:02d}"
        assay = "WES_Twist_Exome" if idx <= 6 else "WGS_PCR_Free"
        batch = "BATCH_NVX_01" if idx % 2 == 1 else "BATCH_NVX_02"
        sex = "XX" if idx % 2 == 0 else "XY"

        samplesheet_rows.append(
            {
                "sample_id": sid,
                "batch_id": batch,
                "assay_type": assay,
                "sex": sex,
                "fastq_r1": f"data/fastq/{sid}_R1_001.fastq.gz",
                "fastq_r2": f"data/fastq/{sid}_R2_001.fastq.gz",
            }
        )

        # Base read count (Millions of read pairs)
        raw_reads_m = round(rng.uniform(52.0, 68.0), 2) if assay == "WES_Twist_Exome" else round(rng.uniform(340.0, 410.0), 2)

        # Inject controlled anomalies for realistic MultiQC SLA alerting:
        # NGS_S08: WARN (slightly elevated PCR duplication 16.8%, borderline depth)
        # NGS_S11: FAIL (degraded library: low Q30 82.4%, high duplication 24.6%, low coverage 18.4x)
        if sid == "NGS_S11":
            pre_q30 = 76.8
            post_q30 = 82.4
            trimmed_retention = 0.865
            mapped_pct = 91.2
            dup_pct = 24.6
            mean_cov = 18.4
            cov_30x_pct = 41.5
            gc_pct = 53.8
        elif sid == "NGS_S08":
            pre_q30 = 88.2
            post_q30 = 92.6
            trimmed_retention = 0.962
            mapped_pct = 97.8
            dup_pct = 16.8
            mean_cov = 33.5
            cov_30x_pct = 86.2
            gc_pct = 49.4
        else:
            pre_q30 = round(rng.uniform(89.5, 93.8), 2)
            post_q30 = round(min(97.8, pre_q30 + rng.uniform(2.8, 4.2)), 2)
            trimmed_retention = round(rng.uniform(0.968, 0.988), 4)
            mapped_pct = round(rng.uniform(98.4, 99.6), 2)
            dup_pct = round(rng.uniform(5.2, 11.4), 2) if assay == "WES_Twist_Exome" else round(rng.uniform(2.1, 5.8), 2)
            mean_cov = round(rng.uniform(95.0, 138.0), 1) if assay == "WES_Twist_Exome" else round(rng.uniform(36.5, 48.0), 1)
            cov_30x_pct = round(rng.uniform(94.2, 98.8), 2) if assay == "WES_Twist_Exome" else round(rng.uniform(89.5, 96.4), 2)
            gc_pct = round(rng.uniform(47.5, 50.8), 2) if assay == "WES_Twist_Exome" else round(rng.uniform(40.5, 42.2), 2)

        trimmed_reads_m = round(raw_reads_m * trimmed_retention, 2)
        mapped_reads_m = round(trimmed_reads_m * (mapped_pct / 100.0), 2)
        dedup_reads_m = round(mapped_reads_m * (1.0 - dup_pct / 100.0), 2)
        insert_size_bp = int(rng.uniform(285, 365))

        # Variant calling yield (WES vs WGS)
        if assay == "WES_Twist_Exome":
            snv_count = int(rng.uniform(38_500, 44_200) * (0.68 if sid == "NGS_S11" else 1.0))
            indel_count = int(snv_count * rng.uniform(0.095, 0.115))
            titv_ratio = round(rng.uniform(2.72, 2.96), 3)
        else:
            snv_count = int(rng.uniform(3_650_000, 3_980_000) * (0.64 if sid == "NGS_S11" else 1.0))
            indel_count = int(snv_count * rng.uniform(0.145, 0.168))
            titv_ratio = round(rng.uniform(2.02, 2.12), 3) if sid != "NGS_S11" else 1.842

        het_hom_ratio = round(rng.uniform(1.52, 1.74), 3)

        telemetry_rows.append(
            {
                "sample_id": sid,
                "batch_id": batch,
                "assay_type": assay,
                "sex": sex,
                "raw_reads_m": raw_reads_m,
                "trimmed_reads_m": trimmed_reads_m,
                "mapped_reads_m": mapped_reads_m,
                "dedup_reads_m": dedup_reads_m,
                "pre_trim_q30_pct": pre_q30,
                "post_trim_q30_pct": post_q30,
                "adapter_trimmed_pct": round((1.0 - trimmed_retention) * 100.0, 2),
                "mapped_pct": mapped_pct,
                "pcr_duplication_pct": dup_pct,
                "median_insert_bp": insert_size_bp,
                "mean_target_cov_x": mean_cov,
                "pct_target_30x": cov_30x_pct,
                "gc_content_pct": gc_pct,
                "pass_snv_count": snv_count,
                "pass_indel_count": indel_count,
                "titv_ratio": titv_ratio,
                "het_hom_ratio": het_hom_ratio,
            }
        )

        # Simulate Nextflow execution trace for the 4 per-sample processes
        scale = 1.0 if assay == "WES_Twist_Exome" else 2.35
        per_sample_processes = [
            ("FASTP_TRIM", 4, 12.0, round(rng.uniform(2.4, 3.8) * (scale ** 0.35), 2), round(rng.uniform(4.5, 8.2) * scale, 2), round(rng.uniform(86.0, 95.5), 1)),
            ("BWA_MEM2_ALIGN", 8, 24.0, round(rng.uniform(14.8, 20.6) * (scale ** 0.15), 2), round(rng.uniform(18.5, 31.0) * scale, 2), round(rng.uniform(91.0, 98.2), 1)),
            ("SAMTOOLS_MARKDUP_STATS", 4, 12.0, round(rng.uniform(4.8, 7.9) * (scale ** 0.25), 2), round(rng.uniform(7.0, 12.8) * scale, 2), round(rng.uniform(82.0, 92.4), 1)),
            ("BCFTOOLS_CALL_FILTER", 4, 12.0, round(rng.uniform(3.9, 6.8) * (scale ** 0.25), 2), round(rng.uniform(9.5, 16.4) * scale, 2), round(rng.uniform(84.5, 94.0), 1)),
        ]
        for proc_name, cpus, alloc_mem, peak_rss, dur_min, cpu_eff in per_sample_processes:
            attempt = 2 if (proc_name == "BWA_MEM2_ALIGN" and sid == "NGS_S09") else 1
            trace_rows.append(
                {
                    "task_id": task_id,
                    "process": proc_name,
                    "tag": sid,
                    "assay_type": assay,
                    "status": "COMPLETED",
                    "exit_code": 0,
                    "attempt": attempt,
                    "cpus": cpus,
                    "allocated_mem_gb": alloc_mem,
                    "peak_rss_gb": round(min(alloc_mem * 0.96, peak_rss), 2),
                    "realtime_min": dur_min,
                    "cpu_efficiency_pct": cpu_eff,
                    "io_read_gb": round(raw_reads_m * 0.14, 2),
                    "io_write_gb": round(raw_reads_m * 0.09, 2),
                    "container": "ghcr.io/denisulko/reproducible-ngs-nextflow:1.0.0",
                }
            )
            task_id += 1

    # Add final cohort aggregation task
    trace_rows.append(
        {
            "task_id": task_id,
            "process": "MULTIQC_COHORT_AGGREGATE",
            "tag": "cohort_multiqc",
            "assay_type": "COHORT_ALL",
            "status": "COMPLETED",
            "exit_code": 0,
            "attempt": 1,
            "cpus": 2,
            "allocated_mem_gb": 4.0,
            "peak_rss_gb": 1.18,
            "realtime_min": 1.45,
            "cpu_efficiency_pct": 88.5,
            "io_read_gb": 0.42,
            "io_write_gb": 0.08,
            "container": "ghcr.io/denisulko/reproducible-ngs-nextflow:1.0.0",
        }
    )

    ss_df = pd.DataFrame(samplesheet_rows)
    tel_df = pd.DataFrame(telemetry_rows)
    tr_df = pd.DataFrame(trace_rows)

    ss_df.to_csv(output_dir / "samplesheet.csv", index=False)
    tel_df.to_csv(output_dir / "pipeline_stage_telemetry.csv", index=False)
    tr_df.to_csv(output_dir / "nextflow_execution_trace.csv", index=False)
    return ss_df, tel_df, tr_df


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/raw")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)

    ss_df, tel_df, tr_df = generate_nextflow_cohort(args.output_dir, seed=args.seed)
    print(f"Generated samplesheet ({len(ss_df)} samples), stage telemetry ({len(tel_df)} rows), and Nextflow trace ({len(tr_df)} tasks).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
