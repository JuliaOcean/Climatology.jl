using Test, Climatology, Statistics, MITgcm, CairoMakie, Suppressor
import NCDatasets, NetCDF, MeshArrays
import Climatology.Dates: Date, DateTime, Month, Day

ENV["DATADEPS_ALWAYS_ACCEPT"]=true

# Environment variables (both optional):
#
#   SKIP_DOWNLOADS=true   bypass downloads / API calls that require server access.
#                         Tests that need downloaded data are skipped (see test_SST.jl).
#   TEST_GROUP=SST        run one group only. A comma-separated list also works, e.g. "SST,ECCO".
#                         Groups: NCDatasetsExt, TimeSeries, SST, ECCO, SLA, SurfaceFluxes.
#                         Default (unset, empty or "all"): run everything.
#
# Examples:
#   TEST_GROUP=SST SKIP_DOWNLOADS=true julia --project -e 'using Pkg; Pkg.test()'
#   TEST_GROUP=SST julia --project=test test/runtests.jl     (after Pkg.develop of this package into test/)
_SKIP_DOWNLOADS = parse(Bool,get(ENV, "SKIP_DOWNLOADS", "false"))

_groups = filter(!isempty,strip.(split(get(ENV, "TEST_GROUP", "all"),",")))
_run(g) = isempty(_groups) || "all" in _groups || g in _groups
let unknown=setdiff(_groups,["all","NCDatasetsExt","TimeSeries","SST","ECCO","SLA","SurfaceFluxes"])
    isempty(unknown) || error("unknown TEST_GROUP: $(join(unknown,", "))")
end

p=dirname(pathof(Climatology))

_run("NCDatasetsExt") && @testset "NCDatasetsExt" begin
    G=MeshArrays.Grids_simple.GridLoad_lonlatdep(collect(1000:1000:5000),ones(360,180,5))
    M=MeshArrays.Integration.define_sums(grid=G,regions=(10,5))
    H=ones(length(M.names),length(M.depths),3)
    V=MeshArrays.Integration.volumes(M,G)
    Climatology.write_H_to_T(tempname()*".nc",M,G,H,V)
end

_run("TimeSeries") && @testset "Time Series Analysis" begin
    dates = collect((DateTime(2020,1,16):Month(1):DateTime(2024,12,31)));
    t_years = Climatology.datetime_to_years(dates);
    dates_reconstructed = Climatology.years_to_datetime(t_years)

    nt=length(dates); data = randn(nt) .+ sin.(2π .* (1:nt) ./ 12)
    data_smc = simple_monthly_climatology(dates, data)
    data_fit = fit_time_series(dates, data, order_season=3)

    @test isa(data_fit,Vector)
end

# Each file below is wrapped in its own @testset and needs only the setup above.
_run("SST") && include("test_SST.jl")
_run("ECCO") && include("test_ECCO.jl")
_run("SLA") && include("test_SLA.jl")

_run("SurfaceFluxes") && @testset "SurfaceFluxes" begin
    data=ERA5.read_sample()
    sst=ERA5.interpolate_sst(data.OISST.sst,data.tim)
    df=ERA5.surface_balance(data.ERA5,sst)
    @test isa(df,Climatology.DataFrame)

    da=Climatology.SurfaceFluxDiag((plot_type=:default,),df)
    plot(da)
    da=Climatology.SurfaceFluxDiag((plot_type=:surface_balance,),(df=df,tim=data.tim,sst=sst))
    plot(da)
    da=Climatology.SurfaceFluxDiag((plot_type=:Qnet_cumsum,),(df=df,tim=data.tim,sst=sst))
    fi=plot(da)
    @test isa(fi,Figure)
end
