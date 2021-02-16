### A Pluto.jl notebook ###
# v0.11.12

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 7b40e150-ee2d-11ea-3346-4b3c11cc782c
using PlutoUI

# ╔═╡ 178d14b0-f5e1-11ea-2c70-dfbdc1b50afc
using Plots

# ╔═╡ c28bac60-edf7-11ea-10ef-c520d1295abc
function U(E,m,x,grad=false)
	a = (1/m)^2
	b = (1/m)^4
	return grad ? (4*E*(-a*(x) + b*(x)^3)) : (4*E*(-(a/2)*(x)^2 + (b/4)*(x)^4))
end

# ╔═╡ 53c36cc0-ee1d-11ea-2786-edebd6fc0d0d
function trapsolver(E, m, viscosity, mass, radius,temp, modamp = 0.0, modf = 0.0)
	T = 1000
	kB = 1.38e-23 # boltzmann constant

	numofsteps = 1000*T
	dt = 1/1000

	x = zeros(Float64, numofsteps)
	v = zeros(Float64, numofsteps)
	
	gamma = 6*π*viscosity*radius # friction coefficient

	for i in 2:numofsteps
		x[i] = x[i-1] + dt*(v[i-1])
		v[i] = v[i-1] + dt*(1/mass)*(-U(E, m, x[i-1], true) - gamma*v[i-1] + modamp*cos(2*pi*modf*dt*(i-1)) + sqrt(2*gamma*kB*temp)*randn())
	end

	t = LinRange(0, T, numofsteps)

	return (x, t)
end

# ╔═╡ 261d4460-ee2e-11ea-20fa-814333077b4c
@bind temp Slider(1:300)

# ╔═╡ b4a749b0-ee2e-11ea-2399-85ac8dbd8e4a
@bind energy_barrier Slider(0:0.00001:1)

# ╔═╡ c521e2a0-ee2e-11ea-3515-3527ebc62271
@bind position Slider(0:0.00001:0.001)

# ╔═╡ 4dc399e0-f5e2-11ea-21d6-89a4a3726108
@bind timelim Slider(0:1000)

# ╔═╡ f1d277b0-ee2e-11ea-354c-a730cba42eb1
x, t = trapsolver(energy_barrier, position, 0.001,1e-6 ,1e-6,temp) # no modulation

# ╔═╡ 1ad7abd0-f5e1-11ea-04ea-b17ce665b0a7
begin
	plot(t,x, xlim=(0,timelim))
end

# ╔═╡ Cell order:
# ╠═c28bac60-edf7-11ea-10ef-c520d1295abc
# ╠═53c36cc0-ee1d-11ea-2786-edebd6fc0d0d
# ╠═7b40e150-ee2d-11ea-3346-4b3c11cc782c
# ╠═261d4460-ee2e-11ea-20fa-814333077b4c
# ╠═b4a749b0-ee2e-11ea-2399-85ac8dbd8e4a
# ╠═c521e2a0-ee2e-11ea-3515-3527ebc62271
# ╠═4dc399e0-f5e2-11ea-21d6-89a4a3726108
# ╠═f1d277b0-ee2e-11ea-354c-a730cba42eb1
# ╠═178d14b0-f5e1-11ea-2c70-dfbdc1b50afc
# ╠═1ad7abd0-f5e1-11ea-04ea-b17ce665b0a7
