# SGLT2i vs. DPP-4i Target Trial Emulation Study Repository

[![Python 3.13](https://img.shields.io/badge/python-3.13-blue.svg)](https://www.python.org/)
[![Status: Contracts Satisfied](https://img.shields.io/badge/Contracts-ALL_SATISFIED-brightgreen.svg)]()
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)]()

This repository contains the complete, reproducible research pipeline, preregistration protocol, and automated data contracts for an observational comparative effectiveness study:  
**"Comparative Cardiovascular Effectiveness and Safety of Sodium-Glucose Cotransporter-2 (SGLT2) Inhibitors versus Dipeptidyl Peptidase-4 (DPP-4) Inhibitors in Patients with Type 2 Diabetes Mellitus: A Preregistered Emulated Target Trial Cohort Study"**

---

## 1. Repository Layout

```
reproducible_cohort_study/
├── data/
│   ├── raw/
│   │   └── synthetic_raw_cohort.csv      # Synthetic EHR fixture with causal DAG
│   ├── interim/
│   │   ├── eligible_cohort.csv           # Cohort after monotonic exclusion filters
│   │   └── attrition_table.csv           # Detailed attrition waterfall ledger
│   └── processed/
│       ├── matched_analysis_cohort.csv   # 1:1 Nearest-neighbor matched cohort
│       └── weighted_analysis_cohort.csv  # Stabilized IPTW analysis cohort
├── notebooks/
│   └── 01_exploratory_and_pipeline_demo.ipynb # Interactive demonstration notebook
├── reports/
│   ├── preregistration.md                # Formal Preregistration Document (PECO, SAP)
│   ├── summary_report.md                 # Execution report & statistical summary
│   ├── pipeline_contract_results.json    # Machine-readable contract verification
│   └── figures/
│       ├── love_plot.png                 # Standardized mean difference balance plot
│       └── km_survival_curves.png        # Kaplan-Meier cumulative incidence curve
├── src/
│   ├── __init__.py
│   ├── config.py                         # Clinical parameters, paths, seeds, thresholds
│   ├── schemas.py                        # Pydantic data contracts and validation models
│   ├── synthetic_generator.py            # Synthetic EHR generator with true causal effect
│   ├── cohort_definition.py              # Monotonic eligibility filters & attrition ledger
│   ├── statistics.py                     # Propensity matching, IPTW, Cox PH, E-values
│   └── data_pipeline.py                  # End-to-end pipeline coordinator
├── tests/
│   ├── __init__.py
│   └── test_pipeline_contracts.py        # Automated data-blind contract test suite
├── environment.lock.yml                  # Cross-platform / Conda environment lockfile
├── requirements.txt                      # Pinned Python requirements
├── requirements.lock                     # Exact hash/version lockfile
├── pyproject.toml                        # Packaging metadata and pytest configuration
└── README.md                             # Reproduction documentation
```

---

## 2. Study Specification (PECO Framework)

- **Population ($P$):** Adults aged 40–85 with diagnosed Type 2 Diabetes Mellitus, baseline HbA1c $\ge 7.0\%$, continuous enrollment $\ge 365$ days, $\text{eGFR} \ge 30\text{ mL/min/1.73m}^2$, without prior MACE or baseline exposure to target drug classes in $[-365, -1]$ days.
- **Exposure ($E$):** Outpatient prescription/dispensing of an SGLT2 inhibitor (empagliflozin, dapagliflozin, canagliflozin).
- **Comparator ($C$):** Outpatient prescription/dispensing of a DPP-4 inhibitor (sitagliptin, linagliptin, saxagliptin).
- **Primary Outcome ($O$):** 3-point Major Adverse Cardiovascular Events (composite of non-fatal MI, non-fatal stroke, or cardiovascular mortality).
- **Negative Control Outcome:** Non-vertebral accidental trauma/fracture (biologically orthogonal control).
- **Observation Window ($W$):** Baseline lookback $[t_0 - 365, t_0 - 1]$; on-treatment follow-up from $t_0 + 1$ to earliest of outcome, treatment discontinuation + 30 days, or administrative end of study (730 days).

---

## 3. Fast Reproduction Guide

### Prerequisites
- Python 3.10+ (tested on Python 3.13)
- Windows / Linux / macOS

### Option A: Using Virtual Environment and Pip Lockfile
```bash
# 1. Create virtual environment
python -m venv .venv

# 2. Activate environment
# On Windows PowerShell:
.venv\Scripts\Activate.ps1
# On Linux/macOS:
source .venv/bin/activate

# 3. Install from reproducible lockfile
pip install -r requirements.lock
```

### Option B: Using Conda Environment Lockfile
```bash
conda env create -f environment.lock.yml
conda activate sglt2-target-trial-emulation
```

---

## 4. Running Pipeline and Verification

### 1. Run Data-Blind Automated Contract Test Suite
```bash
pytest -v tests/test_pipeline_contracts.py
```
*Executes all 7 contract tests, verifying raw schema compliance, monotonic attrition, covariate balance ($\text{SMD} < 0.10$), Cox regression output contracts, negative control calibration, and artifact generation.*

### 2. Run the End-to-End Analysis Pipeline Directly
```bash
python -m src.data_pipeline
```
*Outputs analysis datasets into `data/`, publication figures into `reports/figures/`, and machine-readable verification contract into `reports/pipeline_contract_results.json`.*

### 3. Run the Demonstration Notebook
```bash
jupyter notebook notebooks/01_exploratory_and_pipeline_demo.ipynb
```

---

## 5. Summary of Verified Contracts

| Contract Domain | Prespecified Threshold / Rule | Result | Verification Status |
| :--- | :--- | :---: | :---: |
| **Raw Ingestion Schema** | Non-null, clinical boundary checks | Validated | PASSED |
| **Attrition Monotonicity** | Steps strictly non-increasing | Validated | PASSED |
| **Covariate Balance** | Post-matching $\text{SMD} < 0.10$ across all 15 covariates | **Max SMD = 0.0484** | PASSED |
| **Primary Effect Contract** | As-Treated 3-Point MACE | **HR = 0.777 [0.611, 0.988]** | PASSED ($p = 0.0397$) |
| **Sensitivity Analysis** | Intention-to-Treat 3-Point MACE | **HR = 0.766 [0.616, 0.953]** | PASSED ($p = 0.0167$) |
| **E-Value Robustness** | Quantifies unmeasured confounding tolerance | **E-Value = 1.90** | PASSED |
| **Proportional Hazards** | Piecewise hazard ratio constancy | **$p = 0.6115$** | PASSED |
| **Negative Control** | 95% Confidence Interval must cross 1.00 | **HR = 1.076 [0.810, 1.429]** | PASSED ($p = 0.6127$) |
