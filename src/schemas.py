"""
Data Contracts and Schema Definitions
Defines Pydantic models and validation rules for raw, interim, and processed cohorts,
as well as expected statistical output contracts.
"""
from typing import Optional, List, Dict
from pydantic import BaseModel, Field, field_validator


class RawRecordSchema(BaseModel):
    """Data contract for individual patient records in raw synthetic ingestion."""
    patient_id: str
    age: float = Field(ge=18.0, le=105.0)
    female: int = Field(ge=0, le=1)
    baseline_enrollment_days: int = Field(ge=0)
    index_drug: str  # 'SGLT2i', 'DPP4i', or 'OTHER_HYPOGLYCEMIC'
    prior_sglt2i_use: int = Field(ge=0, le=1)
    prior_dpp4i_use: int = Field(ge=0, le=1)
    baseline_hba1c: float = Field(ge=4.0, le=18.0)
    baseline_egfr: float = Field(ge=5.0, le=160.0)
    baseline_sbp: float = Field(ge=70.0, le=240.0)
    baseline_ldl: float = Field(ge=20.0, le=350.0)
    baseline_ascvd: int = Field(ge=0, le=1)
    baseline_hypertension: int = Field(ge=0, le=1)
    baseline_heart_failure: int = Field(ge=0, le=1)
    smoking_current: int = Field(ge=0, le=1)
    t2dm_duration_gt5y: int = Field(ge=0, le=1)
    baseline_metformin: int = Field(ge=0, le=1)
    baseline_statin: int = Field(ge=0, le=1)
    baseline_ace_arb: int = Field(ge=0, le=1)
    baseline_beta_blocker: int = Field(ge=0, le=1)
    prior_mace_baseline: int = Field(ge=0, le=1)
    
    # Follow-up & Outcomes
    fu_days_primary: float = Field(ge=1.0, le=730.0)
    mace_event: int = Field(ge=0, le=1)
    fu_days_itt: float = Field(ge=1.0, le=730.0)
    mace_event_itt: int = Field(ge=0, le=1)
    
    # Negative Control Outcome (accidental trauma / fracture)
    fu_days_neg_control: float = Field(ge=1.0, le=730.0)
    neg_control_event: int = Field(ge=0, le=1)
    
    # Treatment Discontinuation (for as-treated censoring)
    discontinuation_day: Optional[float] = None


class AttritionStepContract(BaseModel):
    """Schema for individual steps in the cohort attrition waterfall."""
    step_id: int
    criterion_name: str
    retained_count: int = Field(ge=0)
    excluded_count: int = Field(ge=0)
    remaining_sglt2i: int = Field(ge=0)
    remaining_dpp4i: int = Field(ge=0)


class AttritionTableContract(BaseModel):
    """Enforces monotonic descent and non-negative counts across the attrition cascade."""
    steps: List[AttritionStepContract]

    @field_validator("steps")
    def verify_monotonic_retention(cls, v: List[AttritionStepContract]) -> List[AttritionStepContract]:
        for i in range(1, len(v)):
            if v[i].retained_count > v[i - 1].retained_count:
                raise ValueError(
                    f"Attrition violated monotonic decrease at step {v[i].step_id}: "
                    f"{v[i].retained_count} > {v[i-1].retained_count}"
                )
        return v


class CovariateBalanceContract(BaseModel):
    """Enforces standardized mean difference (SMD) threshold for matched cohorts."""
    covariate_name: str
    smd_pre_match: float
    smd_post_match: float = Field(lt=0.10, description="Post-match SMD must be strictly < 0.10")


class StatisticalResultContract(BaseModel):
    """Output contract for primary causal effect estimation."""
    analysis_name: str
    hazard_ratio: float = Field(gt=0.0, lt=10.0)
    ci_lower: float = Field(gt=0.0)
    ci_upper: float = Field(gt=0.0)
    p_value: float = Field(ge=0.0, le=1.0)
    e_value_estimate: float = Field(ge=1.0)
    e_value_ci_bound: float = Field(ge=1.0)
    schoenfeld_ph_pvalue: float = Field(ge=0.0, le=1.0)
    
    @field_validator("ci_upper")
    def ci_upper_must_exceed_lower(cls, v: float, info) -> float:
        if "ci_lower" in info.data and v < info.data["ci_lower"]:
            raise ValueError(f"ci_upper ({v}) cannot be less than ci_lower ({info.data['ci_lower']})")
        return v


class PipelineContractVerification(BaseModel):
    """Summary contract for the entire data-blind pipeline execution."""
    raw_record_count: int = Field(gt=0)
    raw_schema_valid: bool
    attrition_monotonic: bool
    eligible_cohort_size: int = Field(gt=0)
    matched_cohort_size: int = Field(gt=0)
    max_post_match_smd: float = Field(lt=0.10)
    primary_effect: StatisticalResultContract
    negative_control: StatisticalResultContract
    robustness_checks_passed: bool
    execution_status: str
