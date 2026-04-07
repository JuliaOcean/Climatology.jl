
using GLM
import Statistics: mean
import DataFrames: DataFrame

function basis_functions(tt::Vector,zz::Vector)
	t=collect(tt)
	nt=length(t)
	tt0=(collect(1:nt).-0.5)./12
	c1=cos.(2pi.*tt0)
	c2=cos.(2*2pi.*tt0)
	c3=cos.(3*2pi.*tt0)
	s1=sin.(2pi.*tt0)
	s2=sin.(2*2pi.*tt0)
	s3=sin.(3*2pi.*tt0)
    t,(c1,c2,c3),(s1,s2,s3)
end

"""
    fit_monthly_time_series(tt::Vector,zz::Vector; order=0)


Fit a linear trend + seasonal harmonics to specified order (up to 3).

```
zz_fit=Climatology.fit_monthly_time_series(tt,zz; order=0)
```    
"""
function fit_monthly_time_series(tt::Vector,zz::Vector; order=0)
    t,c,s=basis_functions(tt,zz)
	data = DataFrame(property=zz,time=t,
		 c1=c[1],c2=c[2],c3=c[3],s1=s[1],s2=s[2],s3=s[3]);
	if order==0
		mdl = lm(@formula(property ~ time), data)	
	elseif order==1
		mdl = lm(@formula(property ~ time + c1 + s1), data)	
	elseif order==2
		mdl = lm(@formula(property ~ time + c1 + s1 + c2 + s2), data)	
	elseif order==3
		mdl = lm(@formula(property ~ time + c1 + s1 + c2 + s2 + c3 + s3), data)	
    else
        error("order must be between 0 and 3")
	end
    predict(mdl,DataFrame(time=t,c1=c[1],c2=c[2],c3=c[3],s1=s[1],s2=s[2],s3=s[3]))
end


"""
    simple_monthly_climatology(tt::Vector,zz::Vector)


Estimate seasonal simply by averaging all Januaries, all Februaries, etc separately.

Note : this works best in the absence of net trends.

```
zz_mc=Climatology.simple_monthly_climatology(tt,zz)
```    
"""
function simple_monthly_climatology(tt::Vector,zz::Vector)
	nt=length(tt)
    ny=Int(floor(nt/12))
	z=zz[:]
	for m in 1:12
		zmean=vec(mean(z[m:12:ny*12,:],dims=1))
		[z[t,:]=zmean for t in m:12:nt]
	end
	z
end
