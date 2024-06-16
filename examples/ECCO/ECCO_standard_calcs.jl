### A Pluto.jl notebook ###
# v0.19.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5671b654-cfa1-11ec-05c2-6b103f4ebc5f
begin
	using Pkg; Pkg.activate()
	using Distributed, Climatology, PlutoUI

	pth0=pwd()
	pth=ECCO.standard_analysis_setup(pth0)
	Pkg.activate(pth)

	"Done with packages"
end

# ╔═╡ c7ff5c99-ee99-4cef-a35b-bdf8e88691eb
TableOfContents()

# ╔═╡ d85d789a-ca6d-4591-b969-5534fd472018
md"""# Standard Diagnostic Computations

Here we interactively compute standard diagnostics plotted in [ECCO\_standard\_plots.jl](https://gaelforget.github.io/Climatology.jl/dev/examples/ECCO_standard_plots.html) from ECCO climatology.
"""

# ╔═╡ 99a4f24c-863a-4939-b422-6791f14f6341
md"""## Main Computation"""

# ╔═╡ d11ef2d8-760e-4b64-aab3-b2f5c4e379b4
md"""## Configuration"""

# ╔═╡ 7a917876-45fb-4e46-8229-28165fc5c469
begin
	sol0="r2"
	list0=ECCO_helpers.standard_list_toml("")
	P0=ECCO_helpers.parameters(pth,sol0,list0[1])

	!isdir(joinpath(pth,P0.sol)) ? mkdir(joinpath(pth,P0.sol)) : nothing
	pth_trsp=joinpath(pth,P0.sol,"ECCO_transport_lines")
	!isdir(pth_trsp) ? ECCO_helpers.transport_lines(P0.Γ,pth_trsp) : nothing

	readdir(pth)
end

# ╔═╡ 6244de02-bea4-43f6-8963-3723ed192b11
begin
	diag_select = @bind diag confirm(Select(1:length(list0),default=4))
	md"""select a diagnostic : $(diag_select)

	- calc = $(list0[diag].calc)
	- nam  = $(list0[diag].nam)
	- lev  = $(list0[diag].lev)
	"""
end

# ╔═╡ cad6d821-6d49-4a4f-a57c-0511d770fa0c
begin
	P=ECCO_helpers.parameters(P0,list0[diag])
	!isdir(P.pth_out) ? mkdir(P.pth_out) : nothing
	ECCO_diagnostics.driver(P)
	readdir(P.pth_out)
end

# ╔═╡ 60a040f5-02e6-41e2-813a-bf91e0a9c681
md"""## Appendix"""

# ╔═╡ Cell order:
# ╟─c7ff5c99-ee99-4cef-a35b-bdf8e88691eb
# ╟─d85d789a-ca6d-4591-b969-5534fd472018
# ╟─99a4f24c-863a-4939-b422-6791f14f6341
# ╟─6244de02-bea4-43f6-8963-3723ed192b11
# ╠═cad6d821-6d49-4a4f-a57c-0511d770fa0c
# ╟─d11ef2d8-760e-4b64-aab3-b2f5c4e379b4
# ╠═7a917876-45fb-4e46-8229-28165fc5c469
# ╟─60a040f5-02e6-41e2-813a-bf91e0a9c681
# ╠═5671b654-cfa1-11ec-05c2-6b103f4ebc5f
