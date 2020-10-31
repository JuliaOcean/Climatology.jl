# OceanStateEstimation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gaelforget.github.io/OceanStateEstimation.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaelforget.github.io/OceanStateEstimation.jl/dev)
[![Build Status](https://travis-ci.org/gaelforget/OceanStateEstimation.jl.svg?branch=master)](https://travis-ci.org/gaelforget/OceanStateEstimation.jl)

[![Codecov](https://codecov.io/gh/gaelforget/OceanStateEstimation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gaelforget/OceanStateEstimation.jl)
[![Coveralls](https://coveralls.io/repos/github/gaelforget/OceanStateEstimation.jl/badge.svg?branch=master)](https://coveralls.io/github/gaelforget/OceanStateEstimation.jl?branch=master)

[![DOI](https://zenodo.org/badge/260376633.svg)](https://zenodo.org/badge/latestdoi/260376633)

Use examples:

```
using OceanStateEstimation

pth=dirname(pathof(OceanStateEstimation))
lst=joinpath(pth,"../examples/OCCA_climatology.csv")
get_from_dataverse(lst,"DFOsflux",OCCAclim_path)
```

or 

```
using OceanStateEstimation, MeshArrays

γ=GridSpec("LatLonCap","./")
tmp=get_ecco_files(γ,"ETAN")
```

_This package is in early developement stage when breaking changes can be expected._
