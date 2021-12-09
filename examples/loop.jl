
@everywhere sol0="r3"
for ff=1:12
 if ff==1
  @everywhere calc="MHT"
 elseif ff==2
  @everywhere calc="trsp"
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
 include("ECCO_standard_analysis.jl")
end

