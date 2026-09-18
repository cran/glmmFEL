glmmFEL 1.0.6 simulation appraisal (2026-09-18)

These are static results. Installing, loading, documenting, or checking glmmFEL
does not run the simulation. Start with help("glmmFEL-benchmarks").

Design
------
Three response/link combinations: Bernoulli-logit, Bernoulli-probit, Poisson-log.
Independent x ~ N(0,1) and cluster effects b ~ N(0,0.5). Fixed slope 0.5;
intercept -0.5 for binary and log(2) for Poisson. Four cluster settings:
small = 30 x 4; rich = 30 x 20; unbalanced = 30 with sizes 3,5,10,20,32
repeated six times; many_small = 60 x 4. Two hundred paired datasets per cell.
The scenario numbering is in scenarios.csv; seeds and source hashes are in
provenance.json. Predictor draws precede random-effect draws, then responses.

Methods
-------
FEL_Laplace_EM is glmmFEL's modal EM, not glmer's marginal Laplace optimizer.
FEL_mean and FEL_full add posterior mean and variance-diagonal corrections.
FEL_MSPL/FEL_RSPL are the package's marginal/restricted pseudo-likelihood fits.
lme4_Laplace and lme4_AGQ9 use glmer with nAGQ=1/9, bobyqa, maxfun=100000.
glmmML_Laplace and glmmML_GHQ20 use normal random effects, maxit=400,
epsilon=1e-8; glmmML has no probit implementation, recorded as not applicable.
MASS_PQL uses glmmPQL, niter=100, with nlme residual scale fixed at one and
msMaxIter=200. It forces ML for its working Gaussian fit.

glmmFEL EM uses max_iter=400, tol=1e-6 and 1e-6 for all stage tolerances.
Its PL fits use 100 outer and 100 inner iterations. Other controls are defaults.
MASS convergence means returning before its 100th outer iteration, a proxy
because its API does not expose an outer convergence flag. The other methods
use their convergence diagnostics; lme4 singular optima are counted separately
from optimizer failure. Boundary means estimated random-effect variance < 1e-6.

Files and denominators
----------------------
summary.csv: one row per family/setting/method/parameter. Main bias, empirical
variance, RMSE, coverage, and re_mse condition on convergence. bias_mcse is the
sample SD of estimation errors divided by sqrt(converged). rmse_mcse uses the
delta method on squared error. convergence_mcse is the binomial standard error.
Coverage is an exploratory normal-Wald diagnostic using reported SEs; coverage_n
counts usable SEs. The package's SEs condition on the fitted variance component.
boundary_rate divides boundary returned fits by all applicable datasets;
warning_rate divides fits with any warning by all applicable datasets.
all_finite_* retains returned finite estimates even from nonconverged fits and
must not be interpreted as validated inference. re_mse averages latent-effect
prediction loss; different methods return posterior means or conditional modes.

paired-comparisons.csv: FE_mean/FE_full versus each comparator on their common
converged datasets. Negative mse_difference favors the FEL method. Consult
n_common and mse_difference_mcse. This pairing does not eliminate selection
when failures depend on the data. The many comparisons are exploratory.

failures.csv: counts of each nonconvergence/error reason and non-applicability.
stress-summary.csv: separate small crossed/signed-membership convergence grid.
sensitivity-summary.csv and sensitivity-paired.csv: post-hoc 2,000-iteration
EM sensitivity on the same first 40 seeds in every cell, with unchanged main-run
comparators. sensitivity-convergence.csv compares 400/2,000 on those identical
datasets; sensitivity-provenance.json records the separate source/raw hashes.
This smaller sensitivity does not replace the prespecified 200-replicate study.
package-versions.csv: comparator/dependency versions. R was 4.6.1 on Windows.
provenance.json: raw-data hash, frozen-input hashes, seed formula, runtime,
and verification that packaging follow-ups preserve the benchmark values.
Raw per-fit data and exact frozen scripts are retained by the maintainer;
the package distribution contains only the compact static results.

Timings include eight-worker contention on a 20-logical-processor workstation;
they are not isolated microbenchmarks. The study tests one variance and effect
size, a normal covariate, and correctly specified normal random effects. It
does not establish performance for rare events, grouped binomial trials,
larger random-effect variances, misspecification, or arbitrary crossed models.
