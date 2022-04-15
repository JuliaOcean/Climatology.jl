
using OceanStateEstimation, JLD2, Glob, SharedArrays

p=OceanStateEstimation.ECCOdiags_path
q="ECCOv4r3_analysis/"

list_1=readdir(p*q)

list=[]
[push!(list,glob("*.jld2",p*q*r)...) for r in list_1];

#tmp2=[]
#for r in tmp1
#    tmp3=glob("*.jld2",p*q*r)
#    !isempty(tmp3) ? push!(tmp2,tmp3...) : nothing
#end

for r in list
    #tmp1=JLD2.load(r)
    tmp1=collect(keys(JLD2.load(r)))
    if length(tmp1)==1
        println(basename(r))
        tmp2=JLD2.load_object(r)
        if isa(tmp2,SharedArrays.SharedArray)
            println("-")
            println(typeof(tmp2))
            println("=")

            #mv(r,r*"_backup")
            #JLD2.save_object(r,collect(tmp2))
        end
        #println(typeof(tmp2))
        #println(typeof(collect(tmp2)))
    end
end
