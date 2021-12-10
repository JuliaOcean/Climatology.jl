
@everywhere include("ECCO_pkg_grid_etc.jl")

@everywhere sol0="r2"

for ff=1:12
 @everywhere nam="unknown"
 if ff==1
  @everywhere calc="trsp"
 elseif ff==2
  @everywhere calc="MHT"
 elseif ff==3
  @everywhere nam="SIarea"
  @everywhere calc="zonmean2d"
 elseif ff==4
  @everywhere nam="MXLDEPTH"
  @everywhere calc="zonmean2d"
 elseif ff==5
  @everywhere nam="SSH"
  @everywhere calc="zonmean2d"
 elseif ff==6
  @everywhere nam="THETA"
  @everywhere calc="zonmean"	
 elseif ff==7
  @everywhere nam="THETA"
  @everywhere calc="glo2d"
 elseif ff==8
  @everywhere nam="THETA"
  @everywhere calc="glo3d"
 elseif ff==9
  @everywhere nam="SALT"
  @everywhere calc="zonmean"
 elseif ff==10
  @everywhere nam="SALT"
  @everywhere calc="glo2d"
 elseif ff==11
  @everywhere nam="SALT"
  @everywhere calc="glo3d"
 elseif ff==12
  @everywhere calc="overturn"
 end

 @everywhere pth_in,pth_out,pth_tmp,sol,nt,list_steps=ECCO_path_etc(sol0,calc,nam)
 !isdir(pth_out) ? mkdir(pth_out) : nothing
 !isdir(pth_tmp) ? mkdir(pth_tmp) : nothing

 include("ECCO_standard_analysis.jl")
end

