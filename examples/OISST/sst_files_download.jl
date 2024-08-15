
using Distributed

@everywhere include("sst_module.jl")

@everywhere begin
    using Downloads
    path=joinpath(tempdir(),"demo_OISST")
    !ispath(path) ? mkdir(path) : nothing
    fil,_=sst_files.file_lists(path=path)
    list=sst_files.read_files_list(path=path)
    list=list[end-10:end,:]
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
            try
              Downloads.download(r.url[1:end-3]*"_preliminary.nc",r.fil[1:end-3]*"_preliminary.nc")
            catch
              println("file not found online : "*r.fil[1:end-3])
            end
          end
        end
    end
end

else

println("no more files to process")

end
