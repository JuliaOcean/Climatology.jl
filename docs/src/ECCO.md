
The workflow presented here is as follows.

- set up recommended for running analyses of `ECCO estimates`.
- run one computation loop on the `ECCO monthly` files.

```@docs
ECCO.standard_analysis_setup
```

Here is an example of parameters `P` to compute zonal mean temperatures at level 5.

```@docs
ECCO_helpers.parameters
```

The computation loop, over all months, can then be carried out as follows.

```@docs
ECCO_diagnostics.driver
```

```@autodocs
Modules = [OceanStateEstimation.ECCO_io]
```
