"""
Reproducible Data Pipeline Orchestrator
Coordinates Raw Ingestion -> Interim Cohort Curation -> Processed Matched Cohort ->
Statistical Estimation -> Visualizations -> Data Contract Verification.
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, Any, Tuple

from .config import (
    RAW_DATA_PATH,
    INTERIM_DATA_PATH,
    ATTRITION_TABLE_PATH,
    PROCESSED_DATA_PATH,
    WEIGHTED_DATA_PATH,
    CONTRACT_RESULTS_PATH,
    FIGURES_DIR,
    ALL_COVARIATES,
    SMD_ACCEPTANCE_THRESHOLD,
    DEFAULT_SAMPLE_SIZE,
    RANDOM_SEED,
)
from .synthetic_generator import generate_synthetic_raw_cohort
from .schemas import (
    RawRecordSchema,
    CovariateBalanceContract,
    PipelineContractVerification,
)
from .cohort_definition import apply_cohort_eligibility
from .statistics import (
    estimate_propensity_scores,
    trim_propensity_scores,
    perform_nearest_neighbor_matching,
    calculate_iptw_weights,
    calculate_smd,
    fit_cox_ph_model,
    calculate_absolute_risk_and_nnt,
)


def run_pipeline(
    n_patients: int = DEFAULT_SAMPLE_SIZE,
    seed: int = RANDOM_SEED,
    generate_new_data: bool = True,
) -> Dict[str, Any]:
    """
    Executes the complete data-blind pipeline check and records output contracts.
    """
    print("[1/6] Ingestion & Raw Fixture Generation...")
    if generate_new_data or not RAW_DATA_PATH.exists():
        df_raw = generate_synthetic_raw_cohort(
            n_patients=n_patients,
            random_seed=seed,
            output_path="default",
        )
    else:
        df_raw = pd.read_csv(RAW_DATA_PATH)

    # Validate raw schema contract on sample
    raw_sample = df_raw.head(20).to_dict(orient="records")
    for r in raw_sample:
        RawRecordSchema.model_validate(r)
    raw_schema_valid = True

    print("[2/6] Applying Monotonic Cohort Exclusion Filters...")
    df_eligible, df_attrition = apply_cohort_eligibility(df_raw, save_interim=True)
    
    # Check monotonicity of attrition
    counts = df_attrition["retained_count"].tolist()
    attrition_monotonic = all(counts[i] <= counts[i - 1] for i in range(1, len(counts)))

    print("[3/6] Propensity Score Estimation & Caliper Matching...")
    df_ps = estimate_propensity_scores(df_eligible, seed=seed)
    df_trimmed = trim_propensity_scores(df_ps)

    # Pre-matching balance
    smd_pre = calculate_smd(df_trimmed, treatment_col="treatment")
    smd_pre_dict = dict(zip(smd_pre["covariate"], smd_pre["smd"]))

    # 1:1 Matching
    df_matched, caliper_width = perform_nearest_neighbor_matching(df_trimmed, seed=seed)
    
    # Save processed matched cohort
    PROCESSED_DATA_PATH.parent.mkdir(parents=True, exist_ok=True)
    df_matched.to_csv(PROCESSED_DATA_PATH, index=False)

    # Post-matching balance
    smd_post = calculate_smd(df_matched, treatment_col="treatment")
    smd_post_dict = dict(zip(smd_post["covariate"], smd_post["smd"]))

    # Validate Covariate Balance Contracts
    balance_contracts = []
    for cov in ALL_COVARIATES:
        c = CovariateBalanceContract(
            covariate_name=cov,
            smd_pre_match=round(float(smd_pre_dict[cov]), 4),
            smd_post_match=round(float(smd_post_dict[cov]), 4),
        )
        balance_contracts.append(c)

    max_post_smd = max(c.smd_post_match for c in balance_contracts)
    print(f"    -> Max post-match SMD: {max_post_smd:.4f} (Threshold: < {SMD_ACCEPTANCE_THRESHOLD})")

    # Secondary: IPTW
    df_weighted = calculate_iptw_weights(df_trimmed)
    df_weighted.to_csv(WEIGHTED_DATA_PATH, index=False)
    smd_iptw = calculate_smd(df_weighted, treatment_col="treatment", weights_col="iptw_weight")

    print("[4/6] Fitting Outcome Models & Robustness Estimations...")
    # 1. Primary Analysis: As-Treated 3-Point MACE on Matched Cohort
    primary_res = fit_cox_ph_model(
        df_matched,
        duration_col="fu_days_primary",
        event_col="mace_event",
        treatment_col="treatment",
        cluster_col="pair_id",
        analysis_name="Primary As-Treated 3-Point MACE",
    )

    # 2. Intention-to-Treat (ITT) Sensitivity Analysis
    itt_res = fit_cox_ph_model(
        df_matched,
        duration_col="fu_days_itt",
        event_col="mace_event_itt",
        treatment_col="treatment",
        cluster_col="pair_id",
        analysis_name="Sensitivity ITT 3-Point MACE",
    )

    # 3. Negative Control Outcome (Accidental Trauma / Fractures)
    neg_res = fit_cox_ph_model(
        df_matched,
        duration_col="fu_days_neg_control",
        event_col="neg_control_event",
        treatment_col="treatment",
        cluster_col="pair_id",
        analysis_name="Negative Control (Trauma/Fracture)",
    )

    # 4. Absolute Risk & NNT
    abs_metrics = calculate_absolute_risk_and_nnt(df_matched)

    print("[5/6] Generating Publication-Ready Diagnostic Figures...")
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    generate_love_plot(smd_pre, smd_post, FIGURES_DIR / "love_plot.png")
    generate_km_plot(df_matched, FIGURES_DIR / "km_survival_curves.png")

    print("[6/6] Verifying Data Contracts & Generating Report...")
    # Robustness checks passed criteria:
    # 1) ITT consistent in direction (HR < 1.0)
    # 2) Negative control 95% CI includes 1.0
    # 3) Schoenfeld residual p > 0.05
    robustness_passed = (
        itt_res.hazard_ratio < 1.0
        and (neg_res.ci_lower <= 1.0 <= neg_res.ci_upper)
        and (primary_res.schoenfeld_ph_pvalue >= 0.05)
    )

    verification = PipelineContractVerification(
        raw_record_count=len(df_raw),
        raw_schema_valid=raw_schema_valid,
        attrition_monotonic=attrition_monotonic,
        eligible_cohort_size=len(df_eligible),
        matched_cohort_size=len(df_matched),
        max_post_match_smd=max_post_smd,
        primary_effect=primary_res,
        negative_control=neg_res,
        robustness_checks_passed=robustness_passed,
        execution_status="ALL_CONTRACTS_SATISFIED" if (max_post_smd < 0.10 and robustness_passed) else "CONTRACT_VIOLATION",
    )

    # Serialize contract artifact
    contract_dict = verification.model_dump()
    contract_dict["secondary_analyses"] = {
        "itt_effect": itt_res.model_dump(),
        "absolute_risk": abs_metrics,
        "caliper_width": round(caliper_width, 4),
        "balance_table": [c.model_dump() for c in balance_contracts],
    }

    with open(CONTRACT_RESULTS_PATH, "w", encoding="utf-8") as f:
        json.dump(contract_dict, f, indent=2)

    print(f"Contract output successfully written to {CONTRACT_RESULTS_PATH}")
    return contract_dict


def generate_love_plot(smd_pre: pd.DataFrame, smd_post: pd.DataFrame, output_path) -> None:
    """Generates a Love plot illustrating covariate balance before vs after matching."""
    merged = pd.merge(smd_pre, smd_post, on="covariate", suffixes=("_pre", "_post"))
    merged = merged.sort_values("smd_pre", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    y_pos = np.arange(len(merged))

    ax.plot(merged["smd_pre"], y_pos, "o", color="#d9534f", label="Unmatched Cohort", markersize=7)
    ax.plot(merged["smd_post"], y_pos, "s", color="#0275d8", label="1:1 Matched Cohort", markersize=7)
    
    for i in range(len(merged)):
        ax.hlines(y=y_pos[i], xmin=merged.loc[i, "smd_post"], xmax=merged.loc[i, "smd_pre"], color="#aaaaaa", linestyle="--", alpha=0.6)

    ax.axvline(x=0.10, color="#d9534f", linestyle=":", linewidth=1.5, label="Balance Threshold (SMD = 0.10)")
    ax.axvline(x=0.00, color="#333333", linestyle="-", linewidth=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(merged["covariate"], fontsize=9)
    ax.set_xlabel("Standardized Mean Difference (SMD)", fontsize=11, fontweight="bold")
    ax.set_title("Covariate Balance Before and After 1:1 Propensity Score Matching", fontsize=12, pad=12, fontweight="bold")
    ax.set_xlim(-0.02, max(0.40, merged["smd_pre"].max() + 0.05))
    ax.legend(loc="lower right", frameon=True, facecolor="#ffffff", framealpha=0.9)
    ax.grid(axis="x", linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def generate_km_plot(df_matched: pd.DataFrame, output_path) -> None:
    """Generates Kaplan-Meier cumulative incidence survival curve."""
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    time_grid = np.linspace(0, 730, 200)

    for treat_val, label, color in [(1, "SGLT2 Inhibitor", "#2e7d32"), (0, "DPP-4 Inhibitor", "#c62828")]:
        sub = df_matched[df_matched["treatment"] == treat_val]
        times = sub["fu_days_primary"].values
        events = sub["mace_event"].values
        
        # Simple Nelson-Aalen / KM step estimate
        cum_hazard = []
        n_at_risk = len(times)
        cur_events = 0
        
        # Sort events
        order = np.argsort(times)
        s_times = times[order]
        s_events = events[order]

        surv = [1.0]
        t_pts = [0.0]
        cur_surv = 1.0
        
        unique_times = np.unique(s_times[s_events == 1])
        for t in unique_times:
            d = np.sum((s_times == t) & (s_events == 1))
            n = np.sum(s_times >= t)
            if n > 0:
                cur_surv *= (1.0 - d / n)
            t_pts.append(t)
            surv.append(cur_surv)

        t_pts.append(730.0)
        surv.append(cur_surv)

        # Plot cumulative incidence = 1 - S(t)
        cum_inc = (1.0 - np.array(surv)) * 100.0
        ax.step(t_pts, cum_inc, where="post", label=label, color=color, linewidth=2.0)

    ax.set_xlabel("Days from Index Date (Follow-Up)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Cumulative Incidence of 3-Point MACE (%)", fontsize=11, fontweight="bold")
    ax.set_title("Kaplan-Meier Cumulative Incidence (On-Treatment 1:1 Matched)", fontsize=12, pad=12, fontweight="bold")
    ax.set_xlim(0, 730)
    ax.legend(loc="upper left", frameon=True, facecolor="#ffffff", framealpha=0.9)
    ax.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


if __name__ == "__main__":
    run_pipeline()
