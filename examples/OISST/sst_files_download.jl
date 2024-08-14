
using Distributed

@everywhere include("sst_module.jl")

@everywhere begin
    using Downloads
    #list=sst_files.read_files_list(file="ersst_to_get_file_list.csv",add_ymd=false)
    list=sst_files.read_files_list(file="oisst_to_get_file_list.csv",add_ymd=false)
    n_per_workwer=Int(ceil(length(list.fil)/nworkers()))
end

if !isempty(list.fil)

@sync @distributed for m in 1:nworkers()
    n0=n_per_workwer*(m-1)+1
    n1=min(n_per_workwer*m,length(list.fil))
    println("$(n0),$(n1)")
    for r in eachrow(list[n0:n1,:])
        !isdir(dirname(r.fil)) ? mkdir(dirname(r.fil)) : nothing
        if !isfile(r.fil)
          println(r.fil)
          try
            Downloads.download(r.url,r.fil)
          catch
            Downloads.download(r.url[1:end-3]*"_preliminary.nc",r.fil[1:end-3]*"_preliminary.nc")
          end
        end
    end
end

else

println("no more files to process")

end
