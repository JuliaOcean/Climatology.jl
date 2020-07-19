# OceanStateEstimation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gaelforget.github.io/OceanStateEstimation.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaelforget.github.io/OceanStateEstimation.jl/dev)
[![Build Status](https://travis-ci.org/gaelforget/OceanStateEstimation.jl.svg?branch=master)](https://travis-ci.org/gaelforget/OceanStateEstimation.jl)

[![Codecov](https://codecov.io/gh/gaelforget/OceanStateEstimation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gaelforget/OceanStateEstimation.jl)
[![Coveralls](https://coveralls.io/repos/github/gaelforget/OceanStateEstimation.jl/badge.svg?branch=master)](https://coveralls.io/github/gaelforget/OceanStateEstimation.jl?branch=master)

[![DOI](https://zenodo.org/badge/260376633.svg)](https://zenodo.org/badge/latestdoi/260376633)

_This package is in early developement stage when breaking changes can be expected._

```
using MeshArrays, OceanStateEstimation
get_grid_if_needed()
pth="../examples/GRID_LLC90/"
γ=GridSpec("LatLonCap",pth)
Γ=GridLoad(γ)
tmp=get_ecco_files(γ)
```
