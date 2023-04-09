## TruncDeath
Estimation of the Survivor Average Causal Effect with Unmeasured Confounding

Please run the R codes in order.

## leukemiaPKU.csv
Data, including 1161 units.

## SACE_obs.R
Main function to estimate the survivor average causal effect.

    Inputs
    Z: Treatment, 1 for active treatment and 0 for control
    S: Survival, 1 for survival and 0 for death
    Y: Outcome
    X: Cvaraites
    V: Substitutional variable
    subset: Which subset of data to use; default is NULL (to use all data)
    r0: Parameter to control the Mahalanobis distance, a positive number; default is NULL (automatic choice)
    link: Link function to leverage the substitutional variable, linear or ss or iss; default is linear (linear link)
    sensitivity: Sensitivity parameter to nondifferential substitution, a number between -1 and 1; default is 0 (assumption satisfied)
    
    Outputs
    sace: Estimated survivor average causal effect using augmented inverse probability weighting
    sacereg: Estimated survivor average causal effect using regression
    se: Standard error (not considering uncertainty of fitted models)
    sec: Conservative standard error

## saceobs_simulation.R
Simulation of estimation and sensitivity analysis. Outputs are bias and boxplots.

## sensitivity_ns.R
Sensitivity analysis for nondifferential substitution

## sensitivity_sr.R
Sensitivity analysis for substitution relevance

## sensitivity_ni.R
Sensitivity analysis for non-interaction

## sensitivity_monotonicity.R
Sensitivity analysis for monotonicity

## inference_simulation.R
Comparison of confidence intervals.

## dataanalysis.R
Real data analysis, including estimation and sensitivity analysis.
