
using Distributed
year0=1992; year1=2011

@everywhere begin
    using DataFrames, CSV, Dates, Glob
    !isdefined(Main,:sst_files) ? include("sst_module.jl") : nothing
    year0=$year0
    year1=$year1

    datname="oisst"
    list=sst_files.read_files_list(file="$(datname)_whole_file_list.csv",path="files_csv/",add_ymd=true)
    lon,lat=sst_files.read_lon_lat(list.fil[1])

    sel=findall([(f.year>=year0 && f.year<=year1) for f in eachrow(list)])
    suf="$(year0)_$(year1)_"
    gdf=groupby(list[sel,:],:month)
end

output_path=tempname()
@everywhere output_path=$output_path
println("output path="*output_path)
mkdir(output_path)

n_per_workwer=Int(ceil(12/nworkers()))
n_per_workwer*nworkers()!==12 ? println("need nworkers to divide 12") : nothing

for varname in ("sst","anom")
    @sync @distributed for m in 1:nworkers()
        for mm in 1:n_per_workwer
            month=(m-1)*n_per_workwer+mm
            tmp=sst_files.monthlymean(gdf,month,varname=varname)
            sst_files.to_monthly_file(tmp,month,varname=varname,output_path=output_path)
        end
    end
end

output_file=sst_files.write_climatology(output_path,year0,year1,lon,lat)
