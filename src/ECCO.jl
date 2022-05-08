
module ECCO

using Pkg
import OceanStateEstimation: pkg_pth

"""
    ECCO.standard_analysis_setup(pth0)

Create a run folder where folder `pth0` will be linked. 

Folder `pth0` (of type `String`) should be the path to the user's ECCO data folder.

```
using OceanStateEstimation, Pkg

pth0="nctiles_monthly" #edit path as needed
pth=ECCO.standard_analysis_setup(pth0)
Pkg.activate(pth)
Pkg.instantiate()

include("ECCO_standard_loop.jl")
```
"""
function standard_analysis_setup(pth0)
	println(pth0)
	
	#1. setup run folder and create link to ECCO data folder
	pth=joinpath(tempdir(),"ECCO_diags_dev"); 
	!isdir(pth) ? mkdir(pth) : nothing
	pth1=joinpath(pth,"ECCOv4r2")
	!isdir(pth1) ? mkdir(pth1) : nothing
	link0=joinpath(pth1,"nctiles_monthly")
	!isfile(link0)&& !islink(link0) ? symlink(pth0,link0) : nothing
	
	#2. copy Project.toml to run folder
	tmp0=pkg_pth
	tmp1=joinpath(tmp0,"..","examples","ECCO","ECCO_standard_Project.toml")
	tmp2=joinpath(pth,"Project.toml")
	!isfile(tmp2) ? cp(tmp1,tmp2) : nothing
		
	return pth
end


end
