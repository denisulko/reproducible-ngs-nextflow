# Nextflow DSL2 Pipeline Execution & MultiQC Cohort Report

- **Workflow Engine:** `Nextflow DSL2 (>=23.10.0)`
- **Container Image:** `ghcr.io/denisulko/reproducible-ngs-nextflow:1.0.0`
- **Samples Processed:** 12 (`6 WES_Twist_Exome`, `6 WGS_PCR_Free`)
- **SLA Gate Summary:** `10 PASS` | `1 WARN` | `1 FAIL`
- **Total Compute Consumption:** `104.00 CPU Core-Hours` across `49 DSL2 tasks`

## 1. Sample-Level MultiQC & SLA Gate Table

| Sample ID | Assay | Raw Reads (M) | Post-Trim Q30 (%) | Mapped (%) | PCR Dup (%) | Mean Cov (x) | >=30x (%) | SNVs | Ti/Tv | SLA Status | Flags |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `NGS_S01` | `WES_Twist_Exome` | 62.2 | 92.8% | 99.3% | 9.4% | 133.4x | 94.6% | 39,746 | 2.726 | **`PASS`** | `NONE` |
| `NGS_S02` | `WES_Twist_Exome` | 53.5 | 93.9% | 99.4% | 9.7% | 118.1x | 98.7% | 43,227 | 2.927 | **`PASS`** | `NONE` |
| `NGS_S03` | `WES_Twist_Exome` | 56.3 | 97.2% | 98.6% | 9.7% | 102.0x | 96.0% | 41,674 | 2.922 | **`PASS`** | `NONE` |
| `NGS_S04` | `WES_Twist_Exome` | 56.2 | 94.2% | 99.1% | 10.8% | 112.2x | 95.2% | 39,018 | 2.746 | **`PASS`** | `NONE` |
| `NGS_S05` | `WES_Twist_Exome` | 56.3 | 95.2% | 98.9% | 11.1% | 132.7x | 95.4% | 43,701 | 2.792 | **`PASS`** | `NONE` |
| `NGS_S06` | `WES_Twist_Exome` | 56.9 | 93.8% | 98.5% | 8.2% | 98.0x | 97.7% | 41,209 | 2.784 | **`PASS`** | `NONE` |
| `NGS_S07` | `WGS_PCR_Free` | 363.7 | 95.2% | 98.5% | 4.4% | 39.1x | 95.8% | 3,728,541 | 2.041 | **`PASS`** | `NONE` |
| `NGS_S08` | `WGS_PCR_Free` | 408.9 | 92.6% | 97.8% | 16.8% | 33.5x | 86.2% | 3,782,865 | 2.106 | **`WARN`** | `ELEVATED_PCR_DUP` |
| `NGS_S09` | `WGS_PCR_Free` | 407.8 | 97.5% | 99.0% | 2.9% | 41.1x | 89.9% | 3,737,517 | 2.066 | **`PASS`** | `NONE` |
| `NGS_S10` | `WGS_PCR_Free` | 375.2 | 96.2% | 98.5% | 2.8% | 43.3x | 94.2% | 3,943,794 | 2.079 | **`PASS`** | `NONE` |
| `NGS_S11` | `WGS_PCR_Free` | 345.1 | 82.4% | 91.2% | 24.6% | 18.4x | 41.5% | 2,546,873 | 1.842 | **`FAIL`** | `LOW_Q30;LOW_MAPPING;HIGH_PCR_DUP;INSUFFICIENT_DEPTH;TITV_OUTLIER` |
| `NGS_S12` | `WGS_PCR_Free` | 340.6 | 96.2% | 99.5% | 2.6% | 37.8x | 90.2% | 3,849,593 | 2.040 | **`PASS`** | `NONE` |

## 2. Nextflow DSL2 Process Resource Benchmarks

| DSL2 Process | Tasks | CPUs | Alloc RAM (GB) | Mean Peak RSS (GB) | RAM Util (%) | Mean CPU Eff (%) | Mean Walltime (min) | Total Core-Hours | Retries |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `FASTP_TRIM` | 12 | 4 | 12.0 | 3.68 | 30.7% | 90.7% | 10.40 | 8.32 | 0 |
| `BWA_MEM2_ALIGN` | 12 | 8 | 24.0 | 19.60 | 81.7% | 94.3% | 40.04 | 64.06 | 1 |
| `SAMTOOLS_MARKDUP_STATS` | 12 | 4 | 12.0 | 7.04 | 58.6% | 87.0% | 16.76 | 13.41 | 0 |
| `BCFTOOLS_CALL_FILTER` | 12 | 4 | 12.0 | 5.95 | 49.6% | 89.9% | 22.69 | 18.16 | 0 |
| `MULTIQC_COHORT_AGGREGATE` | 1 | 2 | 4.0 | 1.18 | 29.5% | 88.5% | 1.45 | 0.05 | 0 |
