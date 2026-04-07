
using GLM, Dates
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
    fit_time_series(tt::Vector, zz::Vector; order=0)

General time series fitting that handles both numeric and DateTime vectors.
Dispatches to appropriate method based on input type.

```
dates = collect((DateTime(2020,1,16):Month(1):DateTime(2024,12,31)));
nt=length(dates); data = randn(nt) .+ sin.(2π .* (1:nt) ./ 12)
zz_fit = fit_time_series(dates, data; order=1)
```    
"""
function fit_time_series(tt::Vector{DateTime}, zz::Vector; order=0)
	t_years = datetime_to_years(tt)
	fit_time_series(t_years, zz; order=order)
end

function fit_time_series(tt::Vector, zz::Vector; order=0)
    t, c, s = basis_functions(tt, zz)
	data = DataFrame(property=zz, time=t,
		 c1=c[1], c2=c[2], c3=c[3], s1=s[1], s2=s[2], s3=s[3])
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
    predict(mdl, DataFrame(time=t, c1=c[1], c2=c[2], c3=c[3], s1=s[1], s2=s[2], s3=s[3]))
end

"""
    simple_monthly_climatology(tt::Vector, zz::Vector)

Estimate seasonal component by averaging all Januaries, all Februaries, etc separately.
Accepts either numeric vectors or DateTime vectors.

Note: this works best in the absence of net trends.

```
dates = collect((DateTime(2020,1,16):Month(1):DateTime(2024,12,31)));
nt=length(dates); data = randn(nt) .+ sin.(2π .* (1:nt) ./ 12)
zz_mc = simple_monthly_climatology(dates, data)
```    
"""
function simple_monthly_climatology(tt::Vector{DateTime}, zz::Vector)
	t_years = datetime_to_years(tt)
	simple_monthly_climatology(t_years, zz)
end

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

"""
    datetime_to_years(tt::Vector{DateTime})

Convert DateTime vector to time in years, measured from January 1st of the first year in tt.
Returns time values offset so that year values are preserved (first year y0 is added back).

```
t_DT = collect((DateTime(2020,1,16):Month(1):DateTime(2024,12,31)))
t_years = datetime_to_years(t_DT)
# First value ≈ 2020.0, last value ≈ 2024.92
```
"""
function datetime_to_years(tt::Vector{DateTime})
	y0 = year(tt[1])
	t0 = DateTime(y0, 1, 1)  # January 1st of first year
	t_years = [(Dates.value(t - t0) / (1000 * 86400 * 365.25)) for t in tt]
	y0 .+ t_years
end

"""
    years_to_datetime(t_years::Vector)

Convert years (as returned by datetime_to_years) back to DateTime objects.
Infers the reference year from the first value.

```
dates = collect((DateTime(2020,1,16):Month(1):DateTime(2024,12,31)));
t_years = Climatology.datetime_to_years(dates);
dates_reconstructed = Climatology.years_to_datetime(t_years)
```    
"""
function years_to_datetime(t_years::Vector)
	y0 = Int(floor(t_years[1]))
	t0 = DateTime(y0, 1, 1)
	[t0 + Day(round(Int, (t - y0) * 365.25)) for t in t_years]
end

# Handle Date for fit_time_series
function fit_time_series(tt::Vector{Date}, zz::Vector; order=0)
    fit_time_series(DateTime.(tt), zz; order=order)
end

# Handle Date for simple_monthly_climatology
function simple_monthly_climatology(tt::Vector{Date}, zz::Vector)
    simple_monthly_climatology(DateTime.(tt), zz)
end

# Handle AbstractRange for fit_time_series
function fit_time_series(tt::AbstractRange, zz::Vector; order=0)
	fit_time_series(collect(tt), zz; order=order)
end

# Handle AbstractRange for simple_monthly_climatology
function simple_monthly_climatology(tt::AbstractRange, zz::Vector)
	simple_monthly_climatology(collect(tt), zz)
end

# Handle AbstractRange for years_to_datetime
function years_to_datetime(t_years::AbstractRange)
	years_to_datetime(collect(t_years))
end

