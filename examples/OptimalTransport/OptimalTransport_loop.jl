using Distributed

calc_SatToSat=true
calc_ModToMod=false
calc_ModToSat=false

zm_test_case=true
choice_method="emd2" #only for 2D case

test_methods=false

println(calc_SatToSat)
println(calc_ModToMod)
println(calc_ModToSat)
println(choice_method)
println(zm_test_case)

##

pth_output=joinpath(tempdir(),"OptimalTransport_example")
!isdir(pth_output) ? mkdir(pth_output) : nothing

@everywhere using Distributed, DistributedArrays, SharedArrays
@everywhere using OptimalTransport, Statistics, LinearAlgebra
@everywhere using Tulip, Distances, JLD2, Tables, CSV, DataFrames

#@everywhere Cost=load("examples/example_Cost.jld2")["Cost"]
@everywhere M=Tables.matrix(CSV.read("examples/M.csv",DataFrame))
@everywhere S=Tables.matrix(CSV.read("examples/S.csv",DataFrame))

@everywhere nx=size(M,1)

## functions that use the "zonal sum" test case

@everywhere function ModToMod_MS(i,j)
    Cost=Float64.([abs(i-j) for i in 1:nx, j in 1:nx])
    emd2(M[:,i],M[:,j], Cost, Tulip.Optimizer())
end

@everywhere function SatToSat_MS(i,j)
    Cost=Float64.([abs(i-j) for i in 1:nx, j in 1:nx])
    emd2(S[:,i],S[:,j], Cost, Tulip.Optimizer())
end

@everywhere function ModToSat_MS(i,j)
    Cost=Float64.([abs(i-j) for i in 1:nx, j in 1:nx])
    emd2(M[:,i],S[:,j], Cost, Tulip.Optimizer())
    #ε = 0.01
    #γ = sinkhorn_stabilized_epsscaling(M[:,i],S[:,j], Cost, ε; maxiter=5_000)
    #dot(γ, Cost) #compute optimal cost, directly
end

## functions that use the full 2D case

@everywhere function ModToSat(i,j)
    a=Chl_from_Mod[:,:,i][:]
    b=Chl_from_Sat[:,:,j][:]
    a,b=preprocess_Chl(a,b)

    if choice_method=="sinkhorn2"
        ε = 0.05
        sinkhorn2(a,b, Cost, ε)
    elseif choice_method=="emd2"
        emd2(a,b, Cost, Tulip.Optimizer())
    elseif choice_method=="epsscaling"
        ε = 0.01
        γ = sinkhorn_stabilized_epsscaling(a,b, Cost, ε; maxiter=5_000)
        dot(γ, Cost) #compute optimal cost, directly
    end
end

@everywhere function ModToMod(i,j)
    a=Chl_from_Mod[:,:,i][:]
    b=Chl_from_Mod[:,:,j][:]
    a,b=preprocess_Chl(a,b)

    if choice_method=="sinkhorn2"
        ε = 0.05
        sinkhorn2(a,b, Cost, ε)
    elseif choice_method=="emd2"
        emd2(a,b, Cost, Tulip.Optimizer())
    elseif choice_method=="epsscaling"
        ε = 0.01
        γ = sinkhorn_stabilized_epsscaling(a,b, Cost, ε; maxiter=5_000)
        dot(γ, Cost) #compute optimal cost, directly
    end
end

@everywhere function SatToSat(i,j)
    a=Chl_from_Sat[:,:,i][:]
    b=Chl_from_Sat[:,:,j][:]
    a,b=preprocess_Chl(a,b)

    if choice_method=="sinkhorn2"
        ε = 0.05
        sinkhorn2(a,b, Cost, ε)
    elseif choice_method=="emd2"
        emd2(a,b, Cost, Tulip.Optimizer())
    elseif choice_method=="epsscaling"
        ε = 0.01
        γ = sinkhorn_stabilized_epsscaling(a,b, Cost, ε; maxiter=5_000)
        dot(γ, Cost) #compute optimal cost, directly
    end
end

##

@everywhere include("OptimalTransport_setup.jl")

II=[[i,j] for i in 1:12, j in 1:12][:];

using Random; JJ=shuffle(II);

if calc_ModToMod
    d = SharedArray{Float64}(12,12)
    t0=[time()]
    for kk in 1:36
        @sync @distributed for k in (kk-1)*4 .+ collect(1:4)
         i=JJ[k][1]
         j=JJ[k][2]
         zm_test_case ? d[i,j]=ModToMod_MS(i,j) : d[i,j]=ModToMod(i,j)
        end
        dt=time()-t0[1]
        println("ModToMod $(kk) $(dt)")
        t0[1]=time()
        jldsave(joinpath(pth_output,"ModToMod_$(choice_method).jld2"); d = d.s)
    end
end

if calc_SatToSat
    d = SharedArray{Float64}(12,12)
    t0=[time()]
    for kk in 1:36
        @sync @distributed for k in (kk-1)*4 .+ collect(1:4)
         i=JJ[k][1]
         j=JJ[k][2]
         zm_test_case ? d[i,j]=SatToSat_MS(i,j) : d[i,j]=SatToSat(i,j)
        end
        dt=time()-t0[1]
        println("SatToSat $(kk) $(dt)")
        t0[1]=time()
        jldsave(joinpath(pth_output,"SatToSat.jld2"); d = d.s)
    end
end

if calc_ModToSat
    d = SharedArray{Float64}(12,12)
    t0=[time()]
    for kk in 1:36
        @sync @distributed for k in (kk-1)*4 .+ collect(1:4)
         i=JJ[k][1]
         j=JJ[k][2]
         zm_test_case ? d[i,j]=ModToSat_MS(i,j) : d[i,j]=ModToSat(i,j)
        end
        dt=time()-t0[1]
        println("ModToSat $(kk) $(dt)")
        t0[1]=time()
        jldsave(joinpath(pth_output,"ModToSat.jld2"); d = d.s)
    end
end
    
## function used only for testing several methods at once

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
        jldsave(joinpath(pth_output,"ModToMod_methods.jld2"); d = d.s)
    end
end
