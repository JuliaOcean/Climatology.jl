#STATE/state_3d_set1.0000241020.meta
#    'THETA   ' 'SALT    ' 'DRHODR  '
#TRSP/trsp_3d_set1.0000241020.meta
#    'UVELMASS' 'VVELMASS' 'WVELMASS' 'GM_PsiX ' 'GM_PsiY '
#TRSP/trsp_3d_set3.0000241020.meta
#    'DFxE_TH ' 'DFyE_TH ' 'ADVx_TH ' 'ADVy_TH ' 'DFxE_SLT' 'DFyE_SLT' 'ADVx_SLT' 'ADVy_SLT'

function list_time_steps(pth_in)
    if isdir(joinpath(pth_in,"STATE"))
        list=readdir(joinpath(pth_in,"STATE"))
        list=list[findall(occursin.(Ref("state_3d_set1"),list))]
        list=list[findall(occursin.(Ref("data"),list))]
        [list[i][15:end] for i in 1:length(list)]
    else
        list=[]
    end
    return list
end

function get_nr_nt(pth_in,nam)
    nct_path=joinpath(pth_in,nam)
    lst=readdir(nct_path)
    lst=lst[findall(occursin.(Ref(".nc"),lst))]
    fil1=joinpath(nct_path,lst[1])
    ds=NCTiles.NCDatasets.Dataset(fil1)
    siz=size(ds[nam])
    siz[end-1],siz[end]
end


γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
nr=50
XC=GridLoadVar("XC",γ)
YC=GridLoadVar("YC",γ)
RAC=GridLoadVar("RAC",γ)
hFacC=GridLoadVar("hFacC",γ)
DRF=GridLoadVar("DRF",γ)

nansum(x) = sum(filter(!isnan,x))
nansum(x,y) = mapslices(nansum,x,dims=y)

mskC=hFacC./hFacC

tmp=hFacC./hFacC
tmp=[nansum(tmp[i,j].*RAC[i]) for j in 1:nr, i in eachindex(RAC)]
tot_RAC=nansum(tmp,2)

tmp=[nansum(hFacC[i,j].*RAC[i].*DRF[j]) for j in 1:nr, i in eachindex(RAC)]
tot_VOL=nansum(tmp,2)

## generic read function

function read_monthly(sol,nam,t,list)
    var_list3d=("THETA","SALT","UVELMASS","VVELMASS",
                "ADVx_TH","ADVy_TH","DFxE_TH","DFyE_TH")
    mdsio_list3d=("STATE/state_3d_set1","STATE/state_3d_set1",
        "TRSP/trsp_3d_set1","TRSP/trsp_3d_set1","TRSP/trsp_3d_set2",
        "TRSP/trsp_3d_set2","TRSP/trsp_3d_set2","TRSP/trsp_3d_set2")

    var_list2d=("MXLDEPTH","SIarea","sIceLoad","ETAN")
    mdsio_list2d=("STATE/state_2d_set1","STATE/state_2d_set1",
                  "STATE/state_2d_set1","STATE/state_2d_set1")

    if (sol=="ECCOv4r2_analysis")||(sol=="ECCOv4r3_analysis")
        nct_path=joinpath(pth_in,nam)
        if sum(var_list3d.==nam)==1
            tmp=read_nctiles(nct_path,nam,γ,I=(:,:,:,t))
        else
            tmp=read_nctiles(nct_path,nam,γ,I=(:,:,t))
        end
    elseif (sol=="ECCOv4r4_analysis")
        y0=Int(floor((t-1)/12))+1992
        m0=mod1(t,12)
        nct_path=joinpath(pth_in,nam,string(y0))
        m0<10 ? fil=nam*"_$(y0)_0$(m0).nc" : fil=nam*"_$(y0)_$(m0).nc"
        tmp0=NCTiles.NCDatasets.Dataset(joinpath(nct_path,fil))[nam]
        til0=Tiles(γ,90,90)
        if sum(var_list3d.==nam)==1
            tmp=MeshArray(γ,γ.ioPrec,nr)
            for i in 1:13, k in 1:50
              ff=til0[i].face
              ii=collect(til0[i].i)
              jj=collect(til0[i].j)
              tmp[ff,k][ii,jj]=tmp0[:,:,i,k,1]
            end
            tmp
        else
            tmp=MeshArray(γ,γ.ioPrec)
            for i in 1:13
              ff=til0[i].face
              ii=collect(til0[i].i)
              jj=collect(til0[i].j)
              tmp[ff][ii,jj]=tmp0[:,:,i,1]
            end
            tmp
        end
    else
      if !isempty(findall(var_list3d.==nam))
        ii=findall(var_list3d.==nam)[1];
        fil=mdsio_list3d[ii]
        meta=read_meta(joinpath(pth_in,fil*".0000241020.meta"))
        kk=findall(vec(meta.fldList).==nam)[1]
        tmp=read_mdsio(joinpath(pth_in,fil*list[t][14:end]))[:,:,:,kk]
        tmp=mskC*read(tmp,γ)     
      else
        ii=findall(var_list2d.==nam)[1];
        fil=mdsio_list2d[ii]
        meta=read_meta(joinpath(pth_in,fil*".0000241020.meta"))
        kk=findall(vec(meta.fldList).==nam)[1]
        tmp=read_mdsio(joinpath(pth_in,fil*list[t][14:end]))[:,:,kk]
        tmp=mskC[:,1]*read(tmp,XC)
      end
    end
end

##

#same as LatitudeCircles? split out of LatitudeCircles?
function edge_mask(mskCint::MeshArray)
    mskCint=1.0*mskCint

    #treat the case of blank tiles:
    #mskCint[findall(RAC.==0)].=NaN
    
    mskCplus=exchange(mskCint)

    #edge tracer mask
    mskCedge=similar(mskCint)
    for i in eachindex(mskCedge)
        tmp1=mskCplus[i]
        tmp2=tmp1[2:end-1,1:end-2]+tmp1[2:end-1,3:end]+
            tmp1[1:end-2,2:end-1]+tmp1[3:end,2:end-1]
        mskCedge[i]=1.0*(tmp2.>0).*(tmp1[2:end-1,2:end-1].==0.0)
    end

    #edge velocity mask:
    mskWedge=similar(mskCint)
    mskSedge=similar(mskCint)
    for i in eachindex(mskCedge)
        mskWedge[i]=mskCplus[i][2:end-1,2:end-1] - mskCplus[i][1:end-2,2:end-1]
        mskSedge[i]=mskCplus[i][2:end-1,2:end-1] - mskCplus[i][2:end-1,1:end-2]
    end

    #treat the case of blank tiles:
    #mskCedge[findall(isnan.(mskCedge))].=0.0
    #mskWedge[findall(isnan.(mskWedge))].=0.0
    #mskSedge[findall(isnan.(mskSedge))].=0.0

    return mskCedge,mskWedge,mskSedge
end
