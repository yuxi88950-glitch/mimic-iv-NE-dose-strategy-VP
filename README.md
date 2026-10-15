# Vasopressin initiation: primary analysis

A runnable demonstration of the clone-censor-weight analysis using independently generated **synthetic data**, a reduced covariate set and simplified time terms. No MIMIC-IV records are included; the outputs are illustrative and do not reproduce the manuscript estimates.

At a common time zero, NE is 0.25–<0.35 µg/kg/min. Strategy A initiates VP within 3 hours in this range; strategy B waits until NE reaches 0.35–0.5 µg/kg/min, then initiates within 3 hours in that range. Follow-up is 28 days from time zero.

## Run

Requires R and these packages (tested with R 4.5.1):

```r
install.packages(c("data.table", "survival", "Matrix", "fastglm", "glmnet"))
```

From this directory:

```text
Rscript run.R output 500 4
```

The last two arguments specify bootstrap replicates and workers. Use a new or empty output directory.

`R/data.R` generates data and constructs the clones; `R/models.R` and `R/weights.R` fit time-only numerator and ridge-penalized denominator models and truncate stabilized weights at the 1st/99th percentiles. `R/analysis.R` estimates weighted KM curves, 28-day risk and RMST differences, and a patient-cluster robust Cox HR. Each patient bootstrap refits the weights.

Outputs: `primary.csv`, `bootstrap.csv`, `survival.csv`, `KM.pdf`, weight summaries, and the software environment. All contrasts compare A with B.
