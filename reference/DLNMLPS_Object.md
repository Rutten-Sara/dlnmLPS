# dlnmLPS objects

Objects of class `"dlnmLPS"` as returned by
[`dlnmLPS()`](https://rutten-sara.github.io/dlnmLPS/reference/dlnmLPS.md).

## Components

- model:

  The fitted model.

- xi_mode:

  Posterior mode of the latent field.

- Prec:

  Posterior precision matrix.

- Xpred:

  Model design matrix.

- y:

  Response vector.

- offset:

  Model offset.

- smooth:

  Smooth term information.

- v_mode:

  Posterior mode of log-hyperparameters.

- DIC:

  Deviance Information Criterion.

- `pd DIC`:

  Effective number of parameters for DIC.

- WAIC:

  Widely Applicable Information Criterion.

- `pd WAIC`:

  Effective number of parameters for WAIC.

- CPO:

  Conditional Predictive Ordinates.

- constraints:

  Constraint matrix to prevent identifiability issues.

- types:

  A vector with the type of each model component.
