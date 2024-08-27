module Climatology

pkg_pth=dirname(pathof(Climatology))

#read_Dataset : Placeholder to allow NCDatasets extension, which is activated by `using NCDatasets`.
#read_nctiles_alias : Placeholder to allow MITgcmTools extension, which is activated by `using MITgcmTools`.
function read_Dataset end
function read_nctiles_alias end
function read_mdsio_alias end

#packages that extensions import from Climatology
using Glob, RollingFunctions, JLD2, Statistics
function plot_examples end; export plot_examples

include("types.jl")
include("downloads.jl")
include("ECCO.jl")

import Climatology.downloads: get_ecco_files, get_ecco_variable_if_needed, get_ecco_velocity_if_needed
import Climatology.downloads: get_occa_variable_if_needed, get_occa_velocity_if_needed
import Climatology.downloads: ECCOdiags_add, CBIOMESclim_download, MITPROFclim_download
import DataDeps; import DataDeps: @datadep_str

"""
   examples()

List of examples provided in Climatology.jl (full paths)
"""
function examples()
    nb=joinpath(abspath("/"),split(pathof(Climatology),"/")[2:end-2]...,"examples")
#    ex=glob("*/*.jl",nb)
    ex_known=("CBIOMES_climatology_plot.jl","ECCO_standard_plots.jl",
    "HadIOD_viz.jl","NSLCT_notebook.jl","OptimalTransport_demo.jl")
    ex=[glob("*/"*e,nb)[1] for e in ex_known]
end

export @datadep_str, ECCOdiags_add
export get_ecco_variable_if_needed, get_ecco_velocity_if_needed
export get_occa_variable_if_needed, get_occa_velocity_if_needed

export ECCO, ECCO_helpers, ECCO_io, ECCO_diagnostics, ECCO_procs
export ScratchSpaces

__init__() = begin
    ScratchSpaces.__init__scratch()
    downloads.__init__standard_diags()
end

end # module
