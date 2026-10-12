"""
Cohort Definition & Attrition Cascade Module
Applies prespecified exclusion criteria monotonically and produces an attrition table.
"""
import pandas as pd
from typing import Tuple
from .config import (
    MIN_AGE,
    MAX_AGE,
    MIN_BASELINE_ENROLLMENT_DAYS,
    MIN_EGFR_THRESHOLD,
    EXPOSURE_LABEL,
    COMPARATOR_LABEL,
    INTERIM_DATA_PATH,
    ATTRITION_TABLE_PATH,
)
from .schemas import AttritionTableContract, AttritionStepContract


def apply_cohort_eligibility(
    df_raw: pd.DataFrame,
    save_interim: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Applies preregistered inclusion and exclusion criteria monotonically:
    1. Base Population: Any evaluated patient in raw dataset
    2. Continuous Baseline Enrollment: >= 365 days prior to index date
    3. New-User Washout: No prior dispensing of SGLT2i or DPP4i in baseline window
    4. Target Index Therapy: Initiated either SGLT2i or DPP-4i
    5. Age Eligibility: 40 <= age <= 85
    6. Renal Safety Criterion: baseline eGFR >= 30 mL/min/1.73m²
    7. Incident Outcome Washout: No prior MACE in baseline window
    8. Data Completeness: Valid non-null age and sex records
    
    Returns:
        (df_eligible, df_attrition)
    """
    attrition_steps = []
    current_df = df_raw.copy()

    def record_step(step_id: int, criterion: str, excluded_count: int, df_curr: pd.DataFrame):
        remaining_sglt2i = int((df_curr["index_drug"] == EXPOSURE_LABEL).sum())
        remaining_dpp4i = int((df_curr["index_drug"] == COMPARATOR_LABEL).sum())
        attrition_steps.append(
            AttritionStepContract(
                step_id=step_id,
                criterion_name=criterion,
                retained_count=len(df_curr),
                excluded_count=excluded_count,
                remaining_sglt2i=remaining_sglt2i,
                remaining_dpp4i=remaining_dpp4i,
            )
        )

    # Step 0: Base candidate population
    record_step(0, "Candidate Patients with Antidiabetic Dispensings", 0, current_df)

    # Step 1: Continuous Baseline Enrollment >= 365 days
    mask_enrollment = current_df["baseline_enrollment_days"] >= MIN_BASELINE_ENROLLMENT_DAYS
    excl_enrollment = int((~mask_enrollment).sum())
    current_df = current_df[mask_enrollment].copy()
    record_step(1, f"Continuous baseline enrollment >= {MIN_BASELINE_ENROLLMENT_DAYS} days", excl_enrollment, current_df)

    # Step 2: New-User Washout (no prior SGLT2i or DPP4i)
    mask_washout = (current_df["prior_sglt2i_use"] == 0) & (current_df["prior_dpp4i_use"] == 0)
    excl_washout = int((~mask_washout).sum())
    current_df = current_df[mask_washout].copy()
    record_step(2, "Active new-user washout (no prior SGLT2i/DPP4i in baseline)", excl_washout, current_df)

    # Step 3: Target Class Initiation (SGLT2i or DPP4i)
    mask_target_drug = current_df["index_drug"].isin([EXPOSURE_LABEL, COMPARATOR_LABEL])
    excl_target_drug = int((~mask_target_drug).sum())
    current_df = current_df[mask_target_drug].copy()
    record_step(3, "Initiation of SGLT2i or DPP4i index prescription", excl_target_drug, current_df)

    # Step 4: Age Criteria (40 to 85)
    mask_age = (current_df["age"] >= MIN_AGE) & (current_df["age"] <= MAX_AGE)
    excl_age = int((~mask_age).sum())
    current_df = current_df[mask_age].copy()
    record_step(4, f"Age between {MIN_AGE} and {MAX_AGE} years at index date", excl_age, current_df)

    # Step 5: Renal Safety (eGFR >= 30)
    mask_renal = current_df["baseline_egfr"] >= MIN_EGFR_THRESHOLD
    excl_renal = int((~mask_renal).sum())
    current_df = current_df[mask_renal].copy()
    record_step(5, f"Baseline eGFR >= {MIN_EGFR_THRESHOLD} mL/min/1.73m² (no ESRD/dialysis)", excl_renal, current_df)

    # Step 6: Incident Outcome Washout (No prior MACE in baseline)
    mask_prior_mace = current_df["prior_mace_baseline"] == 0
    excl_prior_mace = int((~mask_prior_mace).sum())
    current_df = current_df[mask_prior_mace].copy()
    record_step(6, "Incident cohort restriction (no prior MACE in baseline)", excl_prior_mace, current_df)

    # Step 7: Data Completeness (non-null age and sex)
    mask_complete = current_df["age"].notna() & current_df["female"].notna()
    excl_complete = int((~mask_complete).sum())
    current_df = current_df[mask_complete].copy()
    record_step(7, "Core demographic data completeness", excl_complete, current_df)

    # Add numeric treatment indicator: 1 = SGLT2i, 0 = DPP4i
    current_df["treatment"] = (current_df["index_drug"] == EXPOSURE_LABEL).astype(int)

    # Validate attrition contract
    attrition_contract = AttritionTableContract(steps=attrition_steps)
    df_attrition = pd.DataFrame([s.model_dump() for s in attrition_contract.steps])

    if save_interim:
        INTERIM_DATA_PATH.parent.mkdir(parents=True, exist_ok=True)
        current_df.to_csv(INTERIM_DATA_PATH, index=False)
        df_attrition.to_csv(ATTRITION_TABLE_PATH, index=False)

    return current_df, df_attrition
