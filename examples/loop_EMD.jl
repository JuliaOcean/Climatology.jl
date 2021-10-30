
calc_SatToSat=true
calc_ModToMod=false

println(calc_SatToSat)
println(calc_ModToMod)

@everywhere using Distributed, DistributedArrays, SharedArrays
@everywhere using OptimalTransport, Statistics, LinearAlgebra

@everywhere function ModToSat(i,j)
    a=Chl_from_Mod[:,:,i][:]
    b=Chl_from_Sat[:,:,j][:]
    a,b=preprocess_Chl(a,b)
    ε = 0.05
    sinkhorn2(a,b, Cost, ε)
end

@everywhere function ModToMod(i,j)
    a=Chl_from_Mod[:,:,i][:]
    b=Chl_from_Mod[:,:,j][:]
    a,b=preprocess_Chl(a,b)
    ε = 0.05
    sinkhorn2(a,b, Cost, ε)
end

@everywhere function SatToSat(i,j)
    a=Chl_from_Mod[:,:,i][:]
    b=Chl_from_Mod[:,:,j][:]
    a,b=preprocess_Chl(a,b)
    ε = 0.05
    sinkhorn2(a,b, Cost, ε)
end

@everywhere include("CBIOMES_climatology_EMD.jl")

II=[[i,j] for i in 1:12, j in 1:12][:];

using Random; JJ=shuffle(II);

#@sync @distributed for k in 
#    i=II[k][1]
#    j=II[k][2]
#    d[i,j]=f1(i,j)
#end

if calc_ModToMod
d = SharedArray{Float64}(12,12)
t0=[time()]
for kk in 1:36
    @sync @distributed for k in (kk-1)*4 .+ collect(1:4)
     i=JJ[k][1]
     j=JJ[k][2]
     d[i,j]=ModToMod(i,j)
    end
    dt=time()-t0[1]
    println("ModToMod $(kk) $(dt)")
    t0[1]=time()
    @save "ModToMod_v0.jld2" d;
end
end

if calc_SatToSat
dd = SharedArray{Float64}(12,12)
t0=[time()]
for kk in 1:36
    @sync @distributed for k in (kk-1)*4 .+ collect(1:4)
     i=JJ[k][1]
     j=JJ[k][2]
     dd[i,j]=SatToSat(i,j)
    end
    dt=time()-t0[1]
    println("SatToSat $(kk) $(dt)")
    t0[1]=time()
    @save "SatToSat_v0.jld2" dd;
end
end

