# Preregistration Protocol: Target Trial Emulation Study

**Title:** Comparative Cardiovascular Effectiveness and Safety of Sodium-Glucose Cotransporter-2 (SGLT2) Inhibitors versus Dipeptidyl Peptidase-4 (DPP-4) Inhibitors in Patients with Type 2 Diabetes Mellitus: A Preregistered Emulated Target Trial Cohort Study  
**Preregistration Protocol ID:** PREREG-EHR-SGLT2-2026-001  
**Protocol Version:** 1.0 (Data-Blind Preregistration Specification)  
**Date of Preregistration:** October 12, 2026  
**Registration Stage:** Formal Study Preregistration prior to unblinding or accessing empirical production data  
**Study Type:** Observational Comparative Effectiveness Research (Target Trial Emulation)  

---

## 1. Study Rationale & Clinical Background

Type 2 Diabetes Mellitus (T2DM) confers a two- to four-fold increased risk of major cardiovascular events. Second-line antidiabetic therapy following metformin typically involves choice between Sodium-Glucose Cotransporter-2 (SGLT2) inhibitors and Dipeptidyl Peptidase-4 (DPP-4) inhibitors. While randomized controlled trials (RCTs) such as EMPA-REG OUTCOME and DECLARE-TIMI 58 demonstrated cardioprotective benefits of SGLT2 inhibitors in enriched, high-risk clinical trial cohorts, real-world comparative effectiveness across broader, routinely managed populations requires rigorous observational target trial emulation to avoid common methodological pitfalls (e.g., immortal time bias, confounding by indication, and healthy user bias).

This study protocol establishes a strictly prespecified, data-blind analysis plan emulating a hypothetical pragmatic target trial comparing active new-users of SGLT2 inhibitors with active new-users of DPP-4 inhibitors.

---

## 2. Target Trial Specification and Emulation Protocol (PECO Framework)

Following the framework of Hernán & Robins (*Annals of Internal Medicine*, 2016), the hypothetical target trial and its observational emulation are specified side-by-side:

| Protocol Component | Target Trial (Hypothetical Pragmatic RCT) | Emulated Cohort Study (Observational Emulation) |
| :--- | :--- | :--- |
| **Eligibility Criteria (Population)** | Adults aged 40–85 with diagnosed T2DM, baseline HbA1c $\ge 7.0\%$, eGFR $\ge 30\text{ mL/min/1.73m}^2$, without active malignancy, end-stage renal disease, or prior MACE within the previous 365 days. | Identical criteria ascertained during a 365-day continuous baseline enrollment lookback period ($[t_0 - 365, t_0 - 1]$). |
| **Intervention / Exposure ($E$)** | Random assignment to initiate an SGLT2 inhibitor (empagliflozin, dapagliflozin, or canagliflozin) at standard outpatient dosage. | Dispensing/prescription of an SGLT2 inhibitor without prior dispensing of any SGLT2i or DPP-4i in the 365-day baseline window (active new-user design). |
| **Active Comparator ($C$)** | Random assignment to initiate a DPP-4 inhibitor (sitagliptin, linagliptin, saxagliptin) at standard outpatient dosage. | Dispensing/prescription of a DPP-4 inhibitor without prior dispensing of any SGLT2i or DPP-4i in the 365-day baseline window. |
| **Assignment Strategy** | Random 1:1 allocation. | Propensity score 1:1 nearest-neighbor caliper matching (caliper = $0.2 \times \text{SD}(\text{logit}(PS))$) and inverse probability of treatment weighting (IPTW). |
| **Observation & Follow-Up Window ($W$)** | Day 1 post-randomization ($t_0 + 1$) until outcome, death, trial termination, or loss to follow-up (up to 730 days / 2 years). | Starts at $t_0 + 1$ (eliminating immortal time bias). Follow-up continues until earliest of: event, non-CV death, treatment discontinuation + 30-day grace period, switching classes, or administrative censoring at 730 days. |
| **Primary Outcome ($O$)** | Time to first 3-point Major Adverse Cardiovascular Event (MACE): composite of non-fatal myocardial infarction, non-fatal ischemic stroke, or cardiovascular death. | Time to first diagnostic code for 3-point MACE in inpatient or emergency encounters during follow-up. |
| **Causal Estimand** | Per-protocol effect (as-treated) and Intention-to-Treat (ITT) effect of SGLT2i vs. DPP-4i on 2-year hazard of 3-point MACE. | Observational analog of per-protocol hazard ratio ($\text{HR}$) and ITT hazard ratio using censored/weighted Cox models. |

---

## 3. Prespecified Hypotheses

* **Null Hypothesis ($H_0$):**
  $$\text{Hazard Ratio (HR)}_{\text{SGLT2i vs DPP-4i for 3-point MACE}} = 1.00$$
  There is no difference in the hazard of incident 3-point MACE between new users of SGLT2 inhibitors and new users of DPP-4 inhibitors.
* **Alternative Hypothesis ($H_1$):**
  $$\text{Hazard Ratio (HR)}_{\text{SGLT2i vs DPP-4i for 3-point MACE}} < 1.00$$
  Initiation of SGLT2 inhibitors is associated with a statistically significant reduction in the hazard of incident 3-point MACE compared with DPP-4 inhibitors.
* **Significance Level:** Two-sided $\alpha = 0.05$. Confidence intervals will be reported at the $95\%$ level ($z_{0.975} = 1.96$).

---

## 4. Exclusion Criteria & Attrition Cascade

Cohort eligibility follows a strictly sequential, monotonic exclusion pipeline. The attrition waterfall records the retained sample size and step attrition count at each stage:

```
[Candidate Patients with Antidiabetic Encounters]
       │
       ▼
[Step 1: Continuous Baseline Enrollment]
  Exclude patients with < 365 days of continuous observation prior to index date
       │
       ▼
[Step 2: Active New-User Washout]
  Exclude patients with any dispensing of SGLT2i or DPP-4i during baseline [-365, -1] days
       │
       ▼
[Step 3: Age Eligibility Criteria]
  Exclude patients aged < 40 years or > 85 years at index date
       │
       ▼
[Step 4: Renal Function / Safety Criteria]
  Exclude baseline eGFR < 30 mL/min/1.73m² or documented dialysis / ESRD
       │
       ▼
[Step 5: Prior MACE Washout (Incident Cohort)]
  Exclude patients with any recorded MACE diagnosis in baseline [-365, -1] days
       │
       ▼
[Step 6: Data Completeness Criteria]
  Exclude patients missing baseline age, sex, or essential clinical vitals
       │
       ▼
[Analysis-Eligible Cohort]
       │
       ├─► SGLT2i Initiators (Exposure Group)
       └─► DPP-4i Initiators (Comparator Group)
```

---

## 5. Covariates, Transformations, & Propensity Score Model

### 5.1 Measured Baseline Covariates
All baseline covariates are assessed strictly within the 365-day pre-index window ($[t_0 - 365, t_0 - 1]$):
1. **Demographics:** Age (continuous years, standardized), Sex (binary: female / male).
2. **Clinical Biomarkers:**
   - Baseline HbA1c ($\%$, continuous, standardized).
   - Baseline eGFR ($\text{mL/min/1.73m}^2$, continuous, standardized).
   - Baseline Systolic Blood Pressure ($\text{mmHg}$, continuous, standardized).
   - Baseline LDL Cholesterol ($\text{mg/dL}$, continuous, standardized).
3. **Comorbidities:**
   - Baseline Atherosclerotic Cardiovascular Disease (ASCVD) history (binary).
   - Baseline Hypertension (binary).
   - Baseline Heart Failure history (binary).
   - Smoking status (binary: current/former vs never).
   - Duration of diagnosed T2DM ($> 5$ years vs $\le 5$ years).
4. **Baseline Concomitant Medications:**
   - Metformin use in baseline (binary).
   - Statin use in baseline (binary).
   - ACE inhibitor / ARB use in baseline (binary).
   - Beta-blocker use in baseline (binary).

### 5.2 Propensity Score Model Specification
A multivariable logistic regression model estimates the conditional probability of receiving SGLT2i versus DPP-4i:
$$e(X) = P(Z = 1 \mid X) = \frac{1}{1 + \exp\left(-\left(\beta_0 + \sum_{k=1}^p \beta_k X_k\right)\right)}$$
where $Z \in \{0, 1\}$ denotes treatment initiation, and $X$ is the vector of prespecified baseline covariates.

### 5.3 Propensity Score Trimming & Matching
1. **Trimming:** Patients in non-overlapping regions of the propensity score distribution (below the 2.5th percentile of the SGLT2i group or above the 97.5th percentile of the DPP-4i group) are excluded to enforce positivity.
2. **1:1 Nearest-Neighbor Caliper Matching:** Greedily matches each SGLT2i initiator to a DPP-4i initiator within a caliper width of:
   $$\text{Caliper} = 0.2 \times \text{SD}\left(\text{logit}(e(X))\right)$$
   without replacement.
3. **Secondary Weighting:** Stabilized Inverse Probability of Treatment Weighting (IPTW) with weight trimming at the 1st and 99th percentiles.

### 5.4 Covariate Balance Criterion
Covariate balance before and after matching is quantified using the Standardized Mean Difference (SMD):
$$\text{SMD} = \frac{|\bar{X}_1 - \bar{X}_0|}{\sqrt{\frac{s_1^2 + s_0^2}{2}}}$$
* **Acceptance Threshold:** An $\text{SMD} < 0.10$ for all measured baseline covariates constitutes successful balance. Any covariate exceeding $0.10$ post-matching triggers multivariable adjustment in the outcome model.

---

## 6. Statistical Analysis Plan & Effect Estimation

### 6.1 Primary Outcome Model
Cox proportional hazards regression models the hazard of 3-point MACE in the matched cohort:
$$h(t \mid Z) = h_0(t) \exp(\beta_{\text{treat}} Z)$$
- **Effect Size:** Point estimate $\widehat{\text{HR}} = \exp(\widehat{\beta}_{\text{treat}})$.
- **Standard Error:** Robust sandwich variance estimator clustering on matched pairs to account for correlation.
- **Hypothesis Test:** Wald test $z = \frac{\widehat{\beta}_{\text{treat}}}{\text{SE}(\widehat{\beta}_{\text{treat}})}$ and corresponding two-tailed $p$-value.

### 6.2 Absolute Effect Measures
- **Cumulative Incidence:** Estimated using the Kaplan-Meier estimator $\widehat{S}(t)$.
- **Absolute Risk Reduction (ARR):**
  $$\text{ARR}(t) = \widehat{R}_{\text{DPP4i}}(t) - \widehat{R}_{\text{SGLT2i}}(t)$$
  evaluated at 12 and 24 months.
- **Number Needed to Treat (NNT):**
  $$\text{NNT}(t) = \frac{1}{\text{ARR}(t)}$$

### 6.3 Proportional Hazards Diagnostics
- Evaluated using Schoenfeld residuals against transformed survival time. A correlation test with $p > 0.05$ confirms the proportionality assumption.

---

## 7. Robustness Checks & Sensitivity Analyses

1. **E-Value for Unmeasured Confounding:**
   Calculated following VanderWeele & Ding (*Annals of Internal Medicine*, 2017):
   $$\text{E-value} = \text{HR} + \sqrt{\text{HR} \times (\text{HR} - 1)}$$
   (for $\text{HR} > 1$, or inverted for $\text{HR} < 1$).
   Evaluates the minimum strength of association that an unmeasured confounder would need with both exposure and outcome to explain away the observed effect.
2. **Negative Control Outcome:**
   - *Outcome:* Non-vertebral accidental fracture or trauma.
   - *Rationale:* Not biologically impacted by glucose-lowering mechanism within 2 years.
   - *Criterion:* Expected $\text{HR}_{\text{neg}} \approx 1.00$ with $95\%$ CI crossing 1.00.
3. **Intention-to-Treat (ITT) Sensitivity Analysis:**
   - Patients are analyzed according to initial prescription regardless of subsequent discontinuation, with follow-up censored only at administrative end of study (730 days) or non-CV death.
4. **Caliper Width Sensitivity:**
   - Re-estimation across alternative calipers: $0.05$, $0.10$, and $0.20 \times \text{SD}(\text{logit}(PS))$.
5. **Prespecified Subgroup Analyses:**
   - Stratified by baseline ASCVD (Established ASCVD vs. Multiple Risk Factors).
   - Stratified by Age ($< 65$ years vs. $\ge 65$ years).

---

## 8. Data-Blind Pipeline Contract

To eliminate researcher degrees of freedom and $p$-hacking, the analysis pipeline is pre-tested and validated on fully synthetic data fixtures that mimic the schemas, constraints, and statistical distributions of the target empirical dataset.

The analysis script outputs a machine-readable data contract validation artifact (`pipeline_contract_results.json`) asserting:
1. `raw_schema_valid`: Raw fixture conforms to Pydantic raw schema.
2. `attrition_monotonic`: Attrition counts decrease strictly monotonically.
3. `balance_achieved`: Maximum post-match SMD $< 0.10$.
4. `primary_effect_valid`: HR $\in (0, 10)$, SE $> 0$, $p \in [0, 1]$.
5. `e_value_valid`: E-value $\ge 1.00$.
6. `schoenfeld_ph_valid`: Schoenfeld residuals test $p > 0.05$.
7. `negative_control_calibrated`: Negative control outcome 95% CI contains 1.00.

**Principal Investigators / Signatures:**
- Lead Epidemiologist & Biostatistician: *Research Team*
- Protocol Version: 1.0
- Date: 2026-10-12
