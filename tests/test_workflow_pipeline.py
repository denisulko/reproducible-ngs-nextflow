"""Unit tests for Nextflow DSL2 Workflow Definition & MultiQC SLA Telemetry Engine."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_multiqc_aggregator import evaluate_sample_qc_gates, summarize_nextflow_resources  # noqa: E402


class TestNextflowWorkflowAndTelemetry(unittest.TestCase):
    """Validate Nextflow DSL2 syntax conventions, SLA quality gates, and trace resource math."""

    def test_nextflow_dsl2_files_and_profiles_exist(self) -> None:
        """Verify main.nf and nextflow.config declare DSL2 processes and container/HPC profiles."""
        main_nf = (ROOT / "main.nf").read_text(encoding="utf-8")
        config_nf = (ROOT / "nextflow.config").read_text(encoding="utf-8")

        self.assertIn("nextflow.enable.dsl = 2", main_nf)
        for proc in [
            "FASTP_TRIM",
            "BWA_MEM2_ALIGN",
            "SAMTOOLS_MARKDUP_STATS",
            "BCFTOOLS_CALL_FILTER",
            "MULTIQC_COHORT_AGGREGATE",
        ]:
            self.assertIn(proc, main_nf)

        for profile in ["docker", "singularity", "conda", "slurm", "awsbatch"]:
            self.assertIn(profile, config_nf)

    def test_sla_gate_identifies_pass_warn_and_fail_samples(self) -> None:
        """Verify clinical SLA gate logic distinguishes PASS, WARN, and FAIL samples."""
        pass_sample = {
            "assay_type": "WES_Twist_Exome",
            "post_trim_q30_pct": 94.5,
            "mapped_pct": 99.1,
            "pcr_duplication_pct": 7.2,
            "mean_target_cov_x": 115.0,
            "pct_target_30x": 96.8,
            "titv_ratio": 2.84,
        }
        warn_sample = {
            "assay_type": "WGS_PCR_Free",
            "post_trim_q30_pct": 92.6,
            "mapped_pct": 97.8,
            "pcr_duplication_pct": 16.8,
            "mean_target_cov_x": 33.5,
            "pct_target_30x": 86.2,
            "titv_ratio": 2.06,
        }
        fail_sample = {
            "assay_type": "WGS_PCR_Free",
            "post_trim_q30_pct": 82.4,
            "mapped_pct": 91.2,
            "pcr_duplication_pct": 24.6,
            "mean_target_cov_x": 18.4,
            "pct_target_30x": 41.5,
            "titv_ratio": 1.84,
        }

        res_pass = evaluate_sample_qc_gates(pass_sample)
        res_warn = evaluate_sample_qc_gates(warn_sample)
        res_fail = evaluate_sample_qc_gates(fail_sample)

        self.assertEqual(res_pass["qc_gate_status"], "PASS")
        self.assertEqual(res_pass["qc_flags"], "NONE")

        self.assertEqual(res_warn["qc_gate_status"], "WARN")
        self.assertIn("ELEVATED_PCR_DUP", str(res_warn["qc_flags"]))

        self.assertEqual(res_fail["qc_gate_status"], "FAIL")
        self.assertIn("LOW_Q30", str(res_fail["qc_flags"]))
        self.assertIn("INSUFFICIENT_DEPTH", str(res_fail["qc_flags"]))

    def test_nextflow_resource_aggregation_calculates_core_hours_and_retries(self) -> None:
        """Verify CPU core-hour calculation and automatic task retry counting."""
        trace = pd.DataFrame(
            [
                {
                    "task_id": 1,
                    "process": "BWA_MEM2_ALIGN",
                    "cpus": 8,
                    "allocated_mem_gb": 24.0,
                    "peak_rss_gb": 18.0,
                    "realtime_min": 30.0,
                    "cpu_efficiency_pct": 94.0,
                    "attempt": 1,
                },
                {
                    "task_id": 2,
                    "process": "BWA_MEM2_ALIGN",
                    "cpus": 8,
                    "allocated_mem_gb": 24.0,
                    "peak_rss_gb": 21.0,
                    "realtime_min": 60.0,
                    "cpu_efficiency_pct": 96.0,
                    "attempt": 2,
                },
            ]
        )
        summary = summarize_nextflow_resources(trace)
        self.assertEqual(len(summary), 1)
        row = summary.iloc[0]
        # (8 * 30 / 60) + (8 * 60 / 60) = 4 + 8 = 12.0 core-hours
        self.assertAlmostEqual(float(row["total_core_hours"]), 12.0, places=2)
        self.assertEqual(int(row["retried_tasks"]), 1)
        self.assertAlmostEqual(float(row["mean_peak_rss_gb"]), 19.5, places=2)


if __name__ == "__main__":
    unittest.main()
