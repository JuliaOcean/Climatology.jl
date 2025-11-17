module ScratchSpaces

using Dataverse, Scratch
using Dataverse.downloads.Downloads

# This will be filled in inside `__init__()`
ECCO = ""
OCCA = ""
CBIOMES = ""
MITprof = ""

# Downloads a resource, stores it within path
function download_dataset(url,path)
    fname = joinpath(path, basename(url))
    if !isfile(fname)
        Downloads.download(url, fname)
    end
    return fname
end

function __init__scratch()
    global ECCO = @get_scratch!("ECCO")
    global OCCA = @get_scratch!("OCCA")
    global CBIOMES = @get_scratch!("CBIOMES")
    global MITprof = @get_scratch!("MITprof")
end

end

##

module downloads

import Climatology: pkg_pth
import Climatology: ScratchSpaces
import Climatology: read_nctiles_alias
using Statistics, MeshArrays
using Dataverse, DataDeps, Glob

## Dataverse Donwloads

"""
    get_ecco_files(γ::gcmgrid,v::String,t=1)

```
using MeshArrays, Climatology, MITgcm
γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
Climatology.get_ecco_variable_if_needed("oceQnet")
tmp=read_nctiles(joinpath(ScratchSpaces.ECCO,"oceQnet/oceQnet"),"oceQnet",γ,I=(:,:,1))
```
"""
function get_ecco_files(γ::gcmgrid,v::String,t=1)
    get_ecco_variable_if_needed(v)
    try
        read_nctiles_alias(joinpath(ScratchSpaces.ECCO,"$v/$v"),"$v",γ,I=(:,:,t))
    catch
        error("failed: call to `read_nctiles`
        This method is provided by `MITgcm.jl`
        and now activated by `using MITgcm` ")
    end
end

"""
    get_ecco_variable_if_needed(v::String)

Download ECCO output for variable `v` to scratch space if needed
"""
function get_ecco_variable_if_needed(v::String)
    lst=Dataverse.file_list("doi:10.7910/DVN/3HPRZI")

    fil=joinpath(ScratchSpaces.ECCO,v,v*".0001.nc")
    if !isfile(fil)
        pth1=joinpath(ScratchSpaces.ECCO,v)
        lst1=findall([v==n[1:end-8] for n in lst.filename])
        !isdir(pth1) ? mkdir(pth1) : nothing
        [Dataverse.file_download(lst,v,pth1) for v in lst.filename[lst1]]
    end
end

"""
    get_ecco_velocity_if_needed()

Download ECCO output for `u,v,w` to scratch space if needed
"""
function get_ecco_velocity_if_needed()
    get_ecco_variable_if_needed("UVELMASS")
    get_ecco_variable_if_needed("VVELMASS")
    get_ecco_variable_if_needed("WVELMASS")
end

"""
    get_occa_variable_if_needed(v::String)

Download OCCA output for variable `v` to scratch space if needed
"""
function get_occa_variable_if_needed(v::String)
    lst=Dataverse.file_list("doi:10.7910/DVN/RNXA2A")
    fil=joinpath(ScratchSpaces.OCCA,v*".0406clim.nc")
    !isfile(fil) ? Dataverse.file_download(lst,v,ScratchSpaces.OCCA) : nothing
end

"""
    get_occa_velocity_if_needed()

Download OCCA output for `u,v,w` to scratch space if needed
"""
function get_occa_velocity_if_needed()
    nams = ("DDuvel","DDvvel","DDwvel","DDtheta","DDsalt")
    [get_occa_variable_if_needed(nam) for nam in nams]
    "done"
end

## zenodo.org and other ownloads

st_d_md(txt="ECCO version 4 release 2") = 
    """
    Dataset: standard analysis for the $(txt) ocean state estimate.
    Authors: Gaël Forget
    """

"""
    unpackDV(filepath)

Like DataDeps's `:unpack` but using `Dataverse.untargz` and remove the `.tar.gz` file.
"""    
function unpackDV(filepath)
    tmp_path=Dataverse.untargz(filepath)
    tmp_path2=joinpath(tmp_path,basename(filepath)[1:end-7])
    tmp_path=(ispath(tmp_path2) ? tmp_path2 : tmp_path)
    if isdir(tmp_path)
        [mv(p,joinpath(dirname(filepath),basename(p))) for p in glob("*",tmp_path)]
        [println(joinpath(dirname(filepath),basename(p))) for p in glob("*",tmp_path)]
        rm(filepath)
    else
        rm(filepath)
        mv(tmp_path,joinpath(dirname(filepath),basename(tmp_path)))
    end
    println("done with unpackDV for "*filepath)
end

"""
    __init__standard_diags()

Register data dependency with DataDep.
"""
function __init__standard_diags()
    register(DataDep("ECCO4R1-stdiags",st_d_md("ECCO4 release 1"),
        ["https://zenodo.org/record/6123262/files/ECCOv4r1_analysis.tar.gz"],
        ["eca4ac04553cb8355b4f26f837d5261b02db23a74934fe2d3d19e7298732881f"],
        post_fetch_method=unpackDV))
    register(DataDep("ECCO4R2-stdiags",st_d_md("ECCO4 release 2"),
        ["https://zenodo.org/record/6123272/files/ECCOv4r2_analysis.tar.gz"],
        ["33338d32863cda25ee43bb2614de252958929693d21a515c04a6fa04a1618bcc"],
        post_fetch_method=unpackDV))
    register(DataDep("ECCO4R3-stdiags",st_d_md("ECCO4 release 3"),
        ["https://zenodo.org/record/6123288/files/ECCOv4r3_analysis.tar.gz"],
        ["313d6f5be1c3634658e6bebf91ae063b1a92b9d4fd3fd2844f9fad25aa1ff56d"],
        post_fetch_method=unpackDV))
    register(DataDep("ECCO4R4-stdiags",st_d_md("ECCO4 release 4"),
        ["https://zenodo.org/record/6123127/files/ECCOv4r4_analysis.tar.gz"],
        ["f4f73ec10e1aceefdef119f111d22e3613b98a4d41e64bf587b795907d5b5011"],
        post_fetch_method=unpackDV))
    register(DataDep("ECCO4R5-stdiags",st_d_md("ECCO4 release 5"),
        ["https://zenodo.org/record/7869067/files/ECCOv4r5_rc2_analysis.tar.gz"],
        ["4df8b59a6464f7a780874dd1ca184ec26adce08faaafb54ccc2bfac37355553d"],
        post_fetch_method=unpackDV))
    register(DataDep("OCCA2HR1-stdiags",st_d_md("OCCA2 historical run 1"),
        ["https://zenodo.org/records/11062685/files/OCCA2HR1_analysis.tar.gz"],
        ["c8470bfb961e50dc061c893fe505a65c16c274840097a4c16a32d78e15038991"],
        post_fetch_method=unpackDV))
    register(DataDep("CBIOMES-clim1","CBIOMES global model climatology",
        ["https://zenodo.org/record/5598417/files/CBIOMES-global-alpha-climatology.nc.tar.gz"],
        ["91b882197cc70809147e0609fe0ea5cdb825c2d5bced187b15b9c13a6c11208e"],
        post_fetch_method=unpackDV))
    register(DataDep("CBIOMES-PML1","CBIOMES global model climatology",
        ["https://rsg.pml.ac.uk/shared_files/brj/CBIOMES_ecoregions/ver_0_2_6/gridded_darwin_montly_clim_360_720_ver_0_2_6.nc"]))
    register(DataDep("MITprof-clim1","MITprof gridded climatologies",
        ["https://zenodo.org/record/5101243/files/gcmfaces_climatologies.tar.gz"],
        ["ae1f6a598f458a8d2b084206140c977b174b0611021283b8cda257596b1bde0c"],
        post_fetch_method=unpackDV))        
    register(DataDep("OISST-stats1","SST climatology and time series",
        ["https://zenodo.org/records/13736355/files/OISST_stats.tar.gz"],
        ["e78850d639396eb52b1fdb8113336ad564090f8e2aeed56ac1d94d6a5e6146de"],
        post_fetch_method=unpackDV))
    register(DataDep("ERA5-OISST","Atmosphere and Ocean state time series",
        ["https://zenodo.org/records/17632724/files/AirSeaHeatFlux.tar.gz"],
        ["5f97f7db4ad7735c69633539d5db22c53a4f0b52863eebcb800f98fc0aacf853"],
        post_fetch_method=unpackDV))        
end

"""
    ECCOdiags_add(nam::String)

Add data to the scratch space folder. Known options for `nam` include 
"release1", "release2", "release3", "release4", "release5", and "OCCA2HR1".

Under the hood this is the same as:

```
using Climatology
datadep"ECCO4R1-stdiags"
datadep"ECCO4R2-stdiags"
datadep"ECCO4R3-stdiags"
datadep"ECCO4R4-stdiags"
datadep"ECCO4R5-stdiags"
datadep"OCCA2HR1-stdiags"
```
"""
function ECCOdiags_add(nam::String)
    withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
        if nam=="release1"||nam=="ECCO4R1"
            datadep"ECCO4R1-stdiags"
        elseif nam=="release2"||nam=="ECCO4R2"
            datadep"ECCO4R2-stdiags"
        elseif nam=="release3"||nam=="ECCO4R3"
            datadep"ECCO4R3-stdiags"
        elseif nam=="release4"||nam=="ECCO4R4"
            datadep"ECCO4R4-stdiags"
        elseif nam=="release5"||nam=="ECCO4R5"
            datadep"ECCO4R5-stdiags"
        elseif nam=="OCCA2HR1"
            datadep"OCCA2HR1-stdiags"
        else
            println("unknown solution")
        end
    end
end

"""
    MITPROFclim_download()

Download lazy artifact to scratch space.
"""   
MITPROFclim_download() = withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
    datadep"MITprof-clim1"
end

"""
    CBIOMESclim_download()

Download lazy artifact to scratch space.
"""
CBIOMESclim_download(nam="clim1") = withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
    if nam=="clim1"
        datadep"CBIOMES-clim1"
    elseif nam=="PML1"
        datadep"CBIOMES-PML1"
    else
        println("unknown data set")
    end
end

"""
    OISSTstats_download()

Download lazy artifact to scratch space.
"""   
OISST_stats_download() = withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
    datadep"OISST-stats1"
end


"""
    ERA5_OISST_download()

Download lazy artifact to scratch space.
"""   
ERA5_OISST_download() = withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
    datadep"ERA5-OISST"
end

end
