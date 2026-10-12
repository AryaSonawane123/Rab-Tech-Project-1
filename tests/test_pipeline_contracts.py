"""
Automated Data-Blind Pipeline Contract Tests
Verifies data contracts, monotonic attrition, covariate balance,
statistical output contracts, and robustness checks.
"""
import json
import pytest
import numpy as np
import pandas as pd
from pathlib import Path

from src.config import (
    ALL_COVARIATES,
    SMD_ACCEPTANCE_THRESHOLD,
    CONTRACT_RESULTS_PATH,
    FIGURES_DIR,
)
from src.synthetic_generator import generate_synthetic_raw_cohort
from src.schemas import (
    RawRecordSchema,
    AttritionTableContract,
    PipelineContractVerification,
)
from src.cohort_definition import apply_cohort_eligibility
from src.statistics import (
    estimate_propensity_scores,
    trim_propensity_scores,
    perform_nearest_neighbor_matching,
    calculate_smd,
    fit_cox_ph_model,
    calculate_e_value,
    calculate_absolute_risk_and_nnt,
)
from src.data_pipeline import run_pipeline


@pytest.fixture(scope="module")
def synthetic_raw_data() -> pd.DataFrame:
    """Fixture providing synthetic raw cohort records."""
    return generate_synthetic_raw_cohort(n_patients=4000, random_seed=123)


def test_raw_schema_conformance(synthetic_raw_data: pd.DataFrame):
    """Asserts all records in raw synthetic fixture conform to Pydantic raw schema."""
    sample = synthetic_raw_data.head(50).to_dict(orient="records")
    for row in sample:
        validated = RawRecordSchema.model_validate(row)
        assert 18.0 <= validated.age <= 105.0
        assert validated.baseline_enrollment_days >= 0
        assert validated.baseline_hba1c >= 4.0
        assert validated.baseline_egfr >= 5.0
        assert validated.mace_event in (0, 1)


def test_attrition_monotonicity_and_positive_retention(synthetic_raw_data: pd.DataFrame):
    """Asserts attrition waterfall decreases strictly monotonically with non-empty cohorts."""
    df_eligible, df_attrition = apply_cohort_eligibility(synthetic_raw_data, save_interim=False)
    
    # Verify non-empty eligible cohort
    assert len(df_eligible) > 0
    assert (df_eligible["treatment"] == 1).sum() > 0
    assert (df_eligible["treatment"] == 0).sum() > 0

    # Verify monotonic decrease in step retained counts
    counts = df_attrition["retained_count"].tolist()
    for i in range(1, len(counts)):
        assert counts[i] <= counts[i - 1], f"Step {i} retention ({counts[i]}) exceeded step {i-1} ({counts[i-1]})"

    # Schema validation
    AttritionTableContract.model_validate({"steps": df_attrition.to_dict(orient="records")})


def test_propensity_score_overlap_and_trimming(synthetic_raw_data: pd.DataFrame):
    """Asserts propensity scores fall in (0, 1) and trimming preserves common support."""
    df_eligible, _ = apply_cohort_eligibility(synthetic_raw_data, save_interim=False)
    df_ps = estimate_propensity_scores(df_eligible, seed=123)

    assert "propensity_score" in df_ps.columns
    assert "logit_ps" in df_ps.columns
    assert (df_ps["propensity_score"] > 0.0).all()
    assert (df_ps["propensity_score"] < 1.0).all()

    df_trimmed = trim_propensity_scores(df_ps)
    assert len(df_trimmed) <= len(df_ps)
    assert len(df_trimmed) > 0.85 * len(df_ps)


def test_covariate_balance_contract(synthetic_raw_data: pd.DataFrame):
    """Asserts 1:1 caliper matching achieves SMD < 0.10 for all baseline covariates."""
    df_eligible, _ = apply_cohort_eligibility(synthetic_raw_data, save_interim=False)
    df_ps = estimate_propensity_scores(df_eligible, seed=123)
    df_trimmed = trim_propensity_scores(df_ps)
    df_matched, caliper = perform_nearest_neighbor_matching(df_trimmed, seed=123)

    # 1:1 matching contract: equal treated and control sample sizes
    n_treated = (df_matched["treatment"] == 1).sum()
    n_control = (df_matched["treatment"] == 0).sum()
    assert n_treated == n_control
    assert n_treated > 500

    # Check SMD contract on matched cohort
    smd_post = calculate_smd(df_matched, treatment_col="treatment")
    max_smd = smd_post["smd"].max()

    print(f"\nMax post-match SMD across {len(smd_post)} covariates: {max_smd:.4f}")
    assert max_smd < SMD_ACCEPTANCE_THRESHOLD, f"Post-match SMD {max_smd} exceeded threshold {SMD_ACCEPTANCE_THRESHOLD}"


def test_cox_proportional_hazards_contract(synthetic_raw_data: pd.DataFrame):
    """Asserts Cox model produces valid effect estimates, SE, CI bounds, and diagnostics."""
    df_eligible, _ = apply_cohort_eligibility(synthetic_raw_data, save_interim=False)
    df_ps = estimate_propensity_scores(df_eligible, seed=123)
    df_trimmed = trim_propensity_scores(df_ps)
    df_matched, _ = perform_nearest_neighbor_matching(df_trimmed, seed=123)

    result = fit_cox_ph_model(
        df_matched,
        duration_col="fu_days_primary",
        event_col="mace_event",
        treatment_col="treatment",
        cluster_col="pair_id",
        analysis_name="Test Cox As-Treated",
    )

    assert 0.0 < result.hazard_ratio < 10.0
    assert 0.0 < result.ci_lower < result.hazard_ratio < result.ci_upper
    assert 0.0 <= result.p_value <= 1.0
    assert result.e_value_estimate >= 1.0
    assert result.schoenfeld_ph_pvalue >= 0.05


def test_negative_control_outcome_contract(synthetic_raw_data: pd.DataFrame):
    """Asserts negative control outcome (accidental trauma) yields null effect covering 1.00."""
    df_eligible, _ = apply_cohort_eligibility(synthetic_raw_data, save_interim=False)
    df_ps = estimate_propensity_scores(df_eligible, seed=123)
    df_trimmed = trim_propensity_scores(df_ps)
    df_matched, _ = perform_nearest_neighbor_matching(df_trimmed, seed=123)

    neg_result = fit_cox_ph_model(
        df_matched,
        duration_col="fu_days_neg_control",
        event_col="neg_control_event",
        treatment_col="treatment",
        cluster_col="pair_id",
        analysis_name="Test Negative Control",
    )

    # 95% Confidence interval must cross 1.0 (null effect expected)
    assert neg_result.ci_lower <= 1.0 <= neg_result.ci_upper, (
        f"Negative control CI [{neg_result.ci_lower}, {neg_result.ci_upper}] failed to cover 1.00"
    )


def test_end_to_end_pipeline_contract_serialization():
    """Runs data-blind pipeline end-to-end, validating serialized contract and figures."""
    contract = run_pipeline(n_patients=5000, seed=42, generate_new_data=True)

    assert CONTRACT_RESULTS_PATH.exists()
    assert (FIGURES_DIR / "love_plot.png").exists()
    assert (FIGURES_DIR / "km_survival_curves.png").exists()

    with open(CONTRACT_RESULTS_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)

    verification = PipelineContractVerification.model_validate(data)
    assert verification.execution_status == "ALL_CONTRACTS_SATISFIED"
    assert verification.max_post_match_smd < SMD_ACCEPTANCE_THRESHOLD
    assert verification.robustness_checks_passed is True
