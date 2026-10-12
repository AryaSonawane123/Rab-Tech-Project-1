"""
Configuration and Study Protocol Definitions
Comparative Effectiveness of SGLT2i vs DPP-4i on 3-point MACE.
"""
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_DATA_PATH = DATA_DIR / "raw" / "synthetic_raw_cohort.csv"
INTERIM_DATA_PATH = DATA_DIR / "interim" / "eligible_cohort.csv"
ATTRITION_TABLE_PATH = DATA_DIR / "interim" / "attrition_table.csv"
PROCESSED_DATA_PATH = DATA_DIR / "processed" / "matched_analysis_cohort.csv"
WEIGHTED_DATA_PATH = DATA_DIR / "processed" / "weighted_analysis_cohort.csv"

REPORTS_DIR = BASE_DIR / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"
CONTRACT_RESULTS_PATH = REPORTS_DIR / "pipeline_contract_results.json"
SUMMARY_REPORT_PATH = REPORTS_DIR / "summary_report.md"

# Random Seeds & Generation Hyperparameters
RANDOM_SEED = 42
DEFAULT_SAMPLE_SIZE = 12000

# PECO Definitions
EXPOSURE_LABEL = "SGLT2i"
COMPARATOR_LABEL = "DPP4i"

# Eligibility Thresholds
MIN_AGE = 40
MAX_AGE = 85
MIN_BASELINE_ENROLLMENT_DAYS = 365
MIN_EGFR_THRESHOLD = 30.0
MAX_FOLLOWUP_DAYS = 730  # 2 years
GRACE_PERIOD_DAYS = 30   # for on-treatment censoring

# Propensity Score & Matching Parameters
CALIPER_SD_MULTIPLIER = 0.20
TRIM_QUANTILE_LOWER = 0.025
TRIM_QUANTILE_UPPER = 0.975
SMD_ACCEPTANCE_THRESHOLD = 0.10

# Baseline Covariates to Balance
COVARIATES_CONTINUOUS = [
    "age",
    "baseline_hba1c",
    "baseline_egfr",
    "baseline_sbp",
    "baseline_ldl",
]

COVARIATES_BINARY = [
    "female",
    "baseline_ascvd",
    "baseline_hypertension",
    "baseline_heart_failure",
    "smoking_current",
    "t2dm_duration_gt5y",
    "baseline_metformin",
    "baseline_statin",
    "baseline_ace_arb",
    "baseline_beta_blocker",
]

ALL_COVARIATES = COVARIATES_CONTINUOUS + COVARIATES_BINARY
