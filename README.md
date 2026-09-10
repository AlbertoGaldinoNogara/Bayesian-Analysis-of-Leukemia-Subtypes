# Bayesian analysis of leukemia subtypes

A hierarchical multivariate normal model, fitted with a hand-written Gibbs sampler,
comparing protein expression across four FAB subtypes of acute myeloid leukemia.

Coursework project — *Bayesian Modelling*, MSc Data Analytics for Business,
Università Cattolica. Joint work with **Tommaso Biganzoli**.

## The question

Reverse-phase protein array measurements give 18 protein markers for each patient.
Patients are labelled by FAB subtype (M0, M1, M2, M4 — the remaining subtypes are
dropped for missingness). Do the subtypes differ in mean protein expression once
the four groups are allowed to borrow strength from each other?

## The model

A multivariate normal likelihood inside a hierarchical normal model, so that each
group mean is itself drawn from a population distribution:

```
y_ij | θ_j, Ω⁻¹  ~  N_p(θ_j, Ω⁻¹)        i = 1..n_j,  j = 1..4
θ_j  | θ_0, T_0⁻¹ ~  N_p(θ_0, T_0⁻¹)
Ω  ~ W_p(a, U)     θ_0 ~ N_p(m_0, K_0)     T_0 ~ W_p(b, W)
```

The joint posterior has no closed form; each full conditional does. That is what
makes Gibbs sampling the natural choice rather than a convenience:

| Parameter | Full conditional |
|---|---|
| Ω | Wishart(a + n, (U + S)⁻¹) |
| θ_0 | Normal, precision K_0 + d·T_0 |
| T_0 | Wishart(b + d, (W + Σ(θ_j − θ_0)(θ_j − θ_0)ᵀ)⁻¹) |
| θ_j | Normal, precision Ω·n_j + T_0 |

Priors are deliberately weakly informative (U = W = 0.01·I, a = b = p, wide
variances), so the data dominate the posterior. 2,000 iterations.

## Findings

**1. The subtypes do not separate.** Comparing the posterior of θ_j across the
four groups for AKT and for BAD, the posterior means overlap substantially. On
this data and this model there is no evidence that mean expression of these
markers differs by FAB subtype.

**2. Shrinkage does most of the visible work.** Group means are pulled towards the
population mean θ_0 — the smaller the group, the harder the pull. The shrinkage
plots show which proteins move and in which direction; a marker below the
shrinkage line has its posterior mean dragged up by the overall mean, and vice
versa. This is the mechanism that makes finding (1) unsurprising rather than
disappointing: with four groups and 18 correlated markers, the hierarchy is doing
what it is designed to do.

**3. The posteriors of Ω and T_0 concentrate near zero**, reflecting the weakly
informative priors — the data, not the prior, are driving the result.

## Repository

```
├── R/gibbs_hierarchical.R     the sampler and the figures, top to bottom
├── data/leukemia.csv          18 protein markers x 256 patients, FAB label in column 1
└── Bayesian Analysis of Leukemia Subtypes.pdf   the written report (model derivation included)
```

## Reproducing

```r
install.packages(c("MASS", "mvtnorm", "ggplot2", "gridExtra", "patchwork"))
source("R/gibbs_hierarchical.R")   # ~2000 Gibbs iterations
```

The script reads `data/leukemia.csv` relative to the repository root, so run it
from there (or open the project in RStudio).

## Data

Semicolon-separated RPPA protein expression values with a FAB subtype label.
Provided as course material; the measurements are de-identified and no patient
attribute other than subtype is present. Rows whose FAB label is blank or outside
{M0, M1, M2, M4} are excluded by the script.

## Known limitations

- Four groups is very few to estimate a between-group covariance T_0 from; the
  hierarchy is only weakly identified.
- Convergence is not formally assessed — no multiple chains, no R-hat, no
  effective sample size. The 2,000 draws are taken at face value.
- The comparison is made marker by marker; with 18 markers examined visually,
  there is no multiplicity control.
- Only M0, M1, M2 and M4 survive the missingness filter, so the conclusion says
  nothing about the subtypes that were dropped.

## Note on the committed code

The script as originally written indexed the group loop with the loop *bound*
rather than the loop *variable*:

```r
for (d in 1:D) {
  bn <- n * omega %*% thetaj[[D]] + ...   # D, not d
  thetaj_post[s, , D] <- thetaJ_post      # D, not d
}
```

Every inner iteration therefore updated group 4 and wrote to slice 4, leaving
slices 1-3 as `NA` -- so the script ran the sampler successfully and then failed
at the first shrinkage plot with `need finite 'ylim' values`. Both indices are now
`d`. With the fix the script runs end to end under R 4.3.3, and

```r
mean(results$thetaj_post[,,1][,4] > results$thetaj_post[,,2][,4])
## 0.5915
```

which is the near-coin-flip the report describes: no meaningful difference in BAD
between M0 and M2. The written conclusions stand; only the code needed correcting.
