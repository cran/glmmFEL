# glmmFEL 1.0.6

* Corrected probit tail derivatives through fourth order, safeguarded posterior
  mode and fixed-effect updates, and kept the E-step distribution fixed during
  the corrected M step. Added independent derivative, tilted-integral,
  quadrature, and Gaussian GLS regression checks.
* Failed inner solves and iteration limits now produce warnings and explicit
  convergence diagnostics. Returned posterior moments and PL working quantities
  are refreshed at the final estimates. Invalid responses, designs, and controls
  receive earlier errors.
  A strictly separating fixed-effect direction in binary data is rejected,
  rather than reporting a small score at an effectively infinite estimate.
* Added exact scalar computation paths for disjoint random-effect designs and
  removed the numerical-Jacobian dependency from estimation. Removed unused
  suggested packages; local benchmark dependencies are documented separately.
* Clarified the modal-EM meaning of `Laplace`, the diagonal scope of `FE_full`
  variance corrections, and the limitations of Hessian standard errors.
  PL fits now return `NA` from `logLik()`; the working Gaussian objective remains
  available as `working_logLik`. Cross-method AIC comparisons are inappropriate.
* Added a static simulation appraisal in `help("glmmFEL-benchmarks")`, covering
  both binary links and Poisson responses, four cluster settings, and comparisons
  with lme4, glmmML, and MASS::glmmPQL. Benchmarks do not run during CRAN checks.

# glmmFEL 1.0.5

* Established a development repository from the supplied 1.0.5 source archive,
  with maintenance checks and a local literature index. No algorithm or public
  API changes were made during repository setup.
