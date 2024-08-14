

using Distributed

@everywhere include("sst_module.jl")

@everywhere begin
    datname="oisst"
    varname="sst"

    list=sst_files.read_files_list(file="$(datname)_whole_file_list.csv",path="files_csv/",add_ymd=false)
    outputfile="lowres_$(datname)_$(varname).csv"; 
    
    ind=coarse_grain.indices(list)

    nt=length(list.fil)
    n_per_workwer=Int(ceil(nt/nworkers()))

    file_root=coarse_grain.file_root("$(varname)_lowres_files","$(varname)_lowres_")
end

isdir(dirname(file_root)) ? mv(dirname(file_root),tempname()) : nothing
mkdir(dirname(file_root))

##

@sync @distributed for m in 1:nworkers()
    n0=n_per_workwer*(m-1)+1
    n1=min(n_per_workwer*m,length(list.fil))
    nnl=coarse_grain.nnl
    println("$(n0),$(n1)")
    for n in n0:n1
        r=list[n,:]
        isfile(r.fil) ? fil=r.fil : fil=r.fil[1:end-3]*"_preliminary.nc"
        if isfile(fil)
            #calculate
            ds=sst_files.Dataset(fil)
            tmp=ds[varname][:,:]
            sst=[coarse_grain.areamean(tmp,ii,jj) for ii in 1:nnl*2, jj in 1:nnl]
            #save to csv
            df=sst_files.DataFrame(i=ind.i,j=ind.j,sst=Float32.(sst[ind.k]))
            tmp=split(basename(r.fil),".")[2]
            sst_files.CSV.write(file_root*tmp*".csv",df)
        end
    end
end

##

coarse_grain.merge_files(outputfile=outputfile,path=dirname(file_root),nam=basename(file_root))
