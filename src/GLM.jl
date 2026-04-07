
using GLM, Dates, StatsModels
import Statistics: mean
import DataFrames: DataFrame

function basis_functions(tt::Vector,zz::Vector)
	t1=collect(tt)
    t2=t1.^2
    t3=t1.^2
	nt=length(t1)
	tt0=(collect(1:nt).-0.5)./12
	c1=cos.(2pi.*tt0)
	c2=cos.(2*2pi.*tt0)
	c3=cos.(3*2pi.*tt0)
	s1=sin.(2pi.*tt0)
	s2=sin.(2*2pi.*tt0)
	s3=sin.(3*2pi.*tt0)
    (t1,t2,t3),(c1,c2,c3),(s1,s2,s3)
end

"""
    build_formula(order_poly::Int, order_season::Int)

Build a GLM formula with polynomial trend and seasonal harmonics.

Arguments:
- order_poly: polynomial degree for time trend (0 to 3)
- order_season: seasonal harmonic order (0 to 3)
"""
function build_formula(order_poly::Int, order_season::Int)
    if order_poly < 0 || order_poly > 3
        error("order_poly must be between 0 and 3")
    end
    if order_season < 0 || order_season > 3
        error("order_season must be between 0 and 3")
    end
    
    # Start with linear time on RHS
    rhs = Term(:time)
    
    # Add higher order polynomial terms
    (order_poly>1) ? (rhs = rhs + Term(:time2)) : nothing
    (order_poly>2) ? (rhs = rhs + Term(:time3)) : nothing
    
    # Add seasonal harmonics
    for i in 1:order_season
        rhs = rhs + Term(Symbol("c$i")) + Term(Symbol("s$i"))
    end
    
    # Create the full formula with LHS ~ RHS
    return Term(:property) ~ rhs
end

"""
    fit_time_series(tt::Vector, zz::Vector; order_poly::Int=1, order_season::Int=1)

General time series fitting that handles both numeric and DateTime vectors.
Dispatches to appropriate method based on input type.

```
dates = collect((DateTime(2020,1,16):Month(1):DateTime(2024,12,31)));
nt=length(dates); data = randn(nt) .+ sin.(2π .* (1:nt) ./ 12)
data_fit = fit_time_series(dates, data; order_season=3)
```    
"""
function fit_time_series(tt::Vector{DateTime}, zz::Vector; order_poly::Int=1, order_season::Int=1)
	t_years = datetime_to_years(tt)
	fit_time_series(t_years, zz; order_poly, order_season)
end

function fit_time_series(tt::Vector, zz::Vector; order_poly::Int=1, order_season::Int=1)
    t, c, s = basis_functions(tt, zz)
	data = DataFrame(property=zz, time=t[1], time2=t[2], time3=t[3],
		 c1=c[1], c2=c[2], c3=c[3], s1=s[1], s2=s[2], s3=s[3])
         
    formula=build_formula(order_poly, order_season)
    mdl = lm(formula, data)	

    predict(mdl, DataFrame(time=t[1], time2=t[2], time3=t[3], 
        c1=c[1], c2=c[2], c3=c[3], s1=s[1], s2=s[2], s3=s[3]))
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

