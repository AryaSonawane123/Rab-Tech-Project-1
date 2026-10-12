"""
Synthetic EHR/Claims Data Generator with Ground-Truth Causal Structure
Generates reproducible cohort fixtures adhering to target trial protocols and schemas.
"""
import numpy as np
import pandas as pd
from typing import Optional
from .config import (
    RANDOM_SEED,
    DEFAULT_SAMPLE_SIZE,
    MIN_AGE,
    MAX_AGE,
    MAX_FOLLOWUP_DAYS,
    GRACE_PERIOD_DAYS,
    RAW_DATA_PATH,
)


def generate_synthetic_raw_cohort(
    n_patients: int = DEFAULT_SAMPLE_SIZE,
    random_seed: int = RANDOM_SEED,
    output_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Generates synthetic raw EHR/claims data simulating type 2 diabetes patients,
    with realistic clinical distributions, confounding by indication, and
    exponential proportional hazards time-to-event outcomes.
    """
    rng = np.random.default_rng(random_seed)

    # 1. Base Demographics
    patient_ids = [f"PT-{i:07d}" for i in range(1, n_patients + 1)]
    age = np.clip(rng.normal(loc=62.5, scale=10.5, size=n_patients), 25.0, 95.0)
    female = rng.binomial(n=1, p=0.46, size=n_patients)

    # Baseline Continuous Enrollment (some will fail the 365-day threshold)
    enrollment_days = rng.choice(
        [rng.integers(30, 364), rng.integers(365, 1200)],
        p=[0.08, 0.92],
        size=n_patients,
    )

    # 2. Baseline Clinical Biomarkers
    # HbA1c (elevated in uncontrolled T2DM, mean ~ 8.2%)
    baseline_hba1c = np.clip(rng.normal(loc=8.2, scale=1.4, size=n_patients), 4.5, 15.0)
    
    # eGFR (mean ~ 76, SD ~ 18, some with CKD stage 4-5 < 30)
    baseline_egfr = np.clip(rng.normal(loc=76.0, scale=19.0, size=n_patients), 10.0, 140.0)
    
    # Blood pressure and lipid profiles
    baseline_sbp = np.clip(rng.normal(loc=136.0, scale=16.0, size=n_patients), 85.0, 210.0)
    baseline_ldl = np.clip(rng.normal(loc=102.0, scale=32.0, size=n_patients), 30.0, 280.0)

    # 3. Baseline Comorbidities & Comedications
    # ASCVD history (higher with older age and SBP)
    ascvd_prob = 1.0 / (1.0 + np.exp(-( -2.2 + 0.04 * (age - 60) + 0.015 * (baseline_sbp - 130) )))
    baseline_ascvd = rng.binomial(1, np.clip(ascvd_prob, 0.05, 0.85))

    # Hypertension
    htn_prob = 1.0 / (1.0 + np.exp(-( -0.5 + 0.03 * (age - 60) + 0.03 * (baseline_sbp - 130) )))
    baseline_hypertension = rng.binomial(1, np.clip(htn_prob, 0.2, 0.95))

    # Heart Failure history
    hf_prob = 1.0 / (1.0 + np.exp(-( -2.8 + 0.05 * (age - 60) + 0.6 * baseline_ascvd - 0.02 * (baseline_egfr - 70) )))
    baseline_heart_failure = rng.binomial(1, np.clip(hf_prob, 0.02, 0.70))

    smoking_current = rng.binomial(1, 0.18, size=n_patients)
    t2dm_duration_gt5y = rng.binomial(1, 0.58, size=n_patients)

    baseline_metformin = rng.binomial(1, 0.82, size=n_patients)
    baseline_statin = rng.binomial(1, np.clip(0.60 + 0.25 * baseline_ascvd, 0.3, 0.95))
    baseline_ace_arb = rng.binomial(1, np.clip(0.55 + 0.20 * baseline_hypertension, 0.2, 0.95))
    baseline_beta_blocker = rng.binomial(1, np.clip(0.20 + 0.35 * baseline_ascvd + 0.25 * baseline_heart_failure, 0.1, 0.9))

    # Prior Washout Violations (small proportion had prior SGLT2i or DPP-4i in lookback)
    prior_sglt2i_use = rng.binomial(1, 0.04, size=n_patients)
    prior_dpp4i_use = rng.binomial(1, 0.05, size=n_patients)

    # Prior MACE during baseline lookback (incident cohort excludes these)
    prior_mace_baseline = rng.binomial(1, np.clip(0.06 * baseline_ascvd, 0.01, 0.30))

    # 4. Treatment Selection (Confounding by Indication)
    # Younger age, higher baseline HbA1c, ASCVD, and HF increase likelihood of SGLT2i vs DPP-4i
    logit_treatment = (
        -0.2
        - 0.02 * (age - 62)
        + 0.35 * (baseline_hba1c - 8.0)
        + 0.50 * baseline_ascvd
        + 0.70 * baseline_heart_failure
        + 0.01 * (baseline_egfr - 75)
        - 0.15 * female
        + 0.25 * baseline_statin
    )
    p_sglt2i = 1.0 / (1.0 + np.exp(-logit_treatment))
    
    # Assign index drug: majority SGLT2i or DPP-4i, small fraction other hypoglycemics
    is_sglt2i = rng.binomial(1, p_sglt2i)
    is_target_class = rng.binomial(1, 0.93, size=n_patients)
    
    index_drug = np.where(
        is_target_class == 0,
        "OTHER_HYPOGLYCEMIC",
        np.where(is_sglt2i == 1, "SGLT2i", "DPP4i")
    )

    # 5. Outcome Generation (Cox Proportional Hazards Data Generating Mechanism)
    # True Causal Effect: SGLT2i reduces 3-point MACE by ~22% (HR = 0.78 => beta_treat = -0.248)
    TRUE_BETA_TREAT = np.log(0.78)

    # Prognostic baseline risk score for MACE
    linear_predictor_mace = (
        0.045 * (age - 60)
        + 0.20 * (baseline_hba1c - 8.0)
        - 0.018 * (baseline_egfr - 75)
        + 0.012 * (baseline_sbp - 130)
        + 0.65 * baseline_ascvd
        + 0.55 * baseline_heart_failure
        + 0.28 * smoking_current
        - 0.20 * baseline_statin
        + TRUE_BETA_TREAT * (index_drug == "SGLT2i").astype(float)
    )

    # Baseline hazard lambda_0 (annual event rate ~ 4.5%)
    LAMBDA_0_MACE = 0.045 / 365.25
    individual_hazard_mace = LAMBDA_0_MACE * np.exp(linear_predictor_mace)
    
    # Generate latent event times via inverse CDF: T = -log(U) / lambda
    u_mace = rng.uniform(0.0001, 0.9999, size=n_patients)
    t_latent_mace = -np.log(u_mace) / individual_hazard_mace

    # Treatment Discontinuation (As-Treated grace period)
    # Fraction discontinuing over 2 years: median adherence ~ 500 days
    has_discontinuation = rng.binomial(1, 0.35, size=n_patients)
    discontinuation_day = np.where(
        has_discontinuation == 1,
        rng.uniform(30.0, 700.0, size=n_patients),
        np.nan,
    )
    
    # Follow-up window censoring:
    # On-treatment primary censoring = min(event_time, discontinuation + 30, administrative 730 days)
    as_treated_censor_time = np.where(
        np.isnan(discontinuation_day),
        MAX_FOLLOWUP_DAYS,
        np.minimum(MAX_FOLLOWUP_DAYS, discontinuation_day + GRACE_PERIOD_DAYS),
    )

    fu_days_primary = np.minimum(t_latent_mace, as_treated_censor_time)
    mace_event = (t_latent_mace <= as_treated_censor_time).astype(int)

    # ITT Follow-up (censoring strictly at administrative 730 days)
    fu_days_itt = np.minimum(t_latent_mace, MAX_FOLLOWUP_DAYS)
    mace_event_itt = (t_latent_mace <= MAX_FOLLOWUP_DAYS).astype(int)

    # 6. Negative Control Outcome (Accidental Trauma / Fractures)
    # Independent of antidiabetic drug mechanism (True beta = 0.0, HR = 1.00)
    linear_predictor_neg = 0.02 * (age - 60) + 0.15 * female + 0.10 * smoking_current
    LAMBDA_0_NEG = 0.035 / 365.25
    individual_hazard_neg = LAMBDA_0_NEG * np.exp(linear_predictor_neg)
    u_neg = rng.uniform(0.0001, 0.9999, size=n_patients)
    t_latent_neg = -np.log(u_neg) / individual_hazard_neg
    fu_days_neg_control = np.minimum(t_latent_neg, MAX_FOLLOWUP_DAYS)
    neg_control_event = (t_latent_neg <= MAX_FOLLOWUP_DAYS).astype(int)

    # Assemble raw DataFrame
    df = pd.DataFrame({
        "patient_id": patient_ids,
        "age": np.round(age, 1),
        "female": female,
        "baseline_enrollment_days": enrollment_days,
        "index_drug": index_drug,
        "prior_sglt2i_use": prior_sglt2i_use,
        "prior_dpp4i_use": prior_dpp4i_use,
        "baseline_hba1c": np.round(baseline_hba1c, 2),
        "baseline_egfr": np.round(baseline_egfr, 1),
        "baseline_sbp": np.round(baseline_sbp, 1),
        "baseline_ldl": np.round(baseline_ldl, 1),
        "baseline_ascvd": baseline_ascvd,
        "baseline_hypertension": baseline_hypertension,
        "baseline_heart_failure": baseline_heart_failure,
        "smoking_current": smoking_current,
        "t2dm_duration_gt5y": t2dm_duration_gt5y,
        "baseline_metformin": baseline_metformin,
        "baseline_statin": baseline_statin,
        "baseline_ace_arb": baseline_ace_arb,
        "baseline_beta_blocker": baseline_beta_blocker,
        "prior_mace_baseline": prior_mace_baseline,
        "fu_days_primary": np.round(fu_days_primary, 1),
        "mace_event": mace_event,
        "fu_days_itt": np.round(fu_days_itt, 1),
        "mace_event_itt": mace_event_itt,
        "fu_days_neg_control": np.round(fu_days_neg_control, 1),
        "neg_control_event": neg_control_event,
        "discontinuation_day": np.round(discontinuation_day, 1),
    })

    if output_path is not None:
        p = RAW_DATA_PATH if output_path == "default" else Path(output_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=False)

    return df


if __name__ == "__main__":
    df_raw = generate_synthetic_raw_cohort(output_path="default")
    print(f"Generated synthetic raw cohort with {len(df_raw)} records at {RAW_DATA_PATH}")
