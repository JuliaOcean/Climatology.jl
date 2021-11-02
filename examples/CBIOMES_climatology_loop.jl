using Distributed

calc_SatToSat=false
calc_ModToMod=false
calc_ModToSat=true
test_methods=false

println(calc_SatToSat)
println(calc_ModToMod)
println(calc_ModToSat)
println(test_methods)

@everywhere using Distributed, DistributedArrays, SharedArrays
@everywhere using OptimalTransport, Statistics, LinearAlgebra
@everywhere using Tulip, Distances, JLD2

@everywhere Cost=load("examples/example_Cost.jld2")["Cost"]

@everywhere function ModToSat(i,j)
    a=Chl_from_Mod[:,:,i][:]
    b=Chl_from_Sat[:,:,j][:]
    a,b=preprocess_Chl(a,b)

    if true #reduce problem size
        a=sum(reshape(a,(120,140)),dims=1)[:]
        b=sum(reshape(b,(120,140)),dims=1)[:]
        Cost=Float64.([abs(i-j) for i in 1:140, j in 1:140])
    end

    #ε = 0.05
    #sinkhorn2(a,b, Cost, ε)
    
    emd2(a,b, Cost, Tulip.Optimizer())

    #ε = 0.01
    #γ = sinkhorn_stabilized_epsscaling(a,b, Cost, ε; maxiter=5_000)
    #dot(γ, Cost) #compute optimal cost, directly
end

@everywhere function ModToMod(i,j)
    a=Chl_from_Mod[:,:,i][:]
    b=Chl_from_Mod[:,:,j][:]
    a,b=preprocess_Chl(a,b)

    if true #reduce problem size
        a=sum(reshape(a,(120,140)),dims=1)[:]
        b=sum(reshape(b,(120,140)),dims=1)[:]
        Cost=Float64.([abs(i-j) for i in 1:140, j in 1:140])
    end

    #ε = 0.05
    #sinkhorn2(a,b, Cost, ε)

    emd2(a,b, Cost, Tulip.Optimizer())

    #ε = 0.01
    #γ = sinkhorn_stabilized_epsscaling(a,b, Cost, ε; maxiter=5_000)
    #dot(γ, Cost) #compute optimal cost, directly
end

@everywhere function SatToSat(i,j)
    a=Chl_from_Sat[:,:,i][:]
    b=Chl_from_Sat[:,:,j][:]
    a,b=preprocess_Chl(a,b)

    if true #reduce problem size
        a=sum(reshape(a,(120,140)),dims=1)[:]
        b=sum(reshape(b,(120,140)),dims=1)[:]
        Cost=Float64.([abs(i-j) for i in 1:140, j in 1:140])
    end

    #ε = 0.05
    #sinkhorn2(a,b, Cost, ε)

    emd2(a,b, Cost, Tulip.Optimizer())

    #ε = 0.01
    #γ = sinkhorn_stabilized_epsscaling(a,b, Cost, ε; maxiter=5_000)
    #dot(γ, Cost) #compute optimal cost, directly
end

@everywhere function ModToMod_methods(i,j,mthd=1)
    a=Chl_from_Mod[:,:,i][:]
    b=Chl_from_Mod[:,:,j][:]
    a,b=preprocess_Chl(a,b)
    
    a=sum(reshape(a,(120,140)),dims=1)[:]
    b=sum(reshape(b,(120,140)),dims=1)[:]
    Cost=Float64.([abs(i-j) for i in 1:140, j in 1:140])

    if mthd==1
        ε = 0.05
        sinkhorn2(a,b, Cost, ε)
    elseif mthd==2
        emd2(a,b, Cost, Tulip.Optimizer())
    elseif mthd==3
        ε = 0.005
        γ = sinkhorn_stabilized(a,b, Cost, ε; maxiter=5_000)
        dot(γ, Cost) #compute optimal cost, directly
    elseif mthd==4
        ε = 0.005
        γ = sinkhorn_stabilized_epsscaling(a,b, Cost, ε; maxiter=5_000)
        dot(γ, Cost) #compute optimal cost, directly
#    elseif mthd==5
#        ε = 0.05
#        γ = quadreg(a,b, Cost, ε; maxiter=100)
#        dot(γ, Cost) #compute optimal cost, directly
end
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
d = SharedArray{Float64}(12,12)
t0=[time()]
for kk in 1:36
    @sync @distributed for k in (kk-1)*4 .+ collect(1:4)
     i=JJ[k][1]
     j=JJ[k][2]
     d[i,j]=SatToSat(i,j)
    end
    dt=time()-t0[1]
    println("SatToSat $(kk) $(dt)")
    t0[1]=time()
    @save "SatToSat_v0.jld2" d;
end
end

if calc_ModToSat
    d = SharedArray{Float64}(12,12)
    t0=[time()]
    for kk in 1:36
        @sync @distributed for k in (kk-1)*4 .+ collect(1:4)
         i=JJ[k][1]
         j=JJ[k][2]
         d[i,j]=ModToSat(i,j)
        end
        dt=time()-t0[1]
        println("ModToSat $(kk) $(dt)")
        t0[1]=time()
        @save "ModToSat_v0.jld2" d;
    end
end
    

if test_methods
    #KK=([1,1],[1,2],[1,9])
    KK=[[1,j] for j in 1:12]
    d = SharedArray{Float64}(6,length(KK))
    t0=[time()]
    for k in 1:4
        for kk in 1:12
         i=KK[kk][1]
         j=KK[kk][2]
         try
            d[k,kk]=ModToMod_methods(i,j,k)
         catch
            d[k,kk]=NaN
         end
         println("$(k) $(kk) $(d[k,kk])")
        end
        dt=time()-t0[1]
        println("ModToMod_methods $(k) $(dt)")
        t0[1]=time()
        @save "ModToMod_methods_v0.jld2" d;
    end
end
