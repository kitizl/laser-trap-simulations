using JLD

function U(a, b, x, grad=false)
	if grad
		return -a*x + b*x^3
	else
		return (-a/2)*x^2 + (b/4)*x^4
	end
end

function trapsolver(a, b, modamp = 0.0, modf = 0.0)
	T = 100000

	numofsteps = 1000*T
	dt = 1/1000

	x = zeros(Float64, numofsteps)
	v = zeros(Float64, numofsteps)

	for i in 2:numofsteps
		x[i] = x[i-1] + dt*(v[i-1])
		v[i] = v[i-1] + dt*(-U(a, b, x[i-1], true) - v[i-1] + modamp*cos(2*pi*modf*dt*(i-1)) + 10*randn())
	end

	t = LinRange(0, T, numofsteps)

	return (t, x)
end

# no modulation
@time t, x = trapsolver(1,16)
save("data-unmod.jld", "t", t, "x", x)
# 1e-3 Hz
@time t, x = trapsolver(1,16,0.01,1e-3)
save("data-mod-1.jld", "t", t, "x", x)
# 1e-2 Hz
@time t, x = trapsolver(1,16,0.01,1e-2)
save("data-mod-2.jld", "t", t, "x", x)
# 1e-1 Hz
@time t, x = trapsolver(1,16,0.01,1e-1)
save("data-mod-3.jld", "t", t, "x", x)
# 1 Hz
@time t, x = trapsolver(1,16,0.01,1.0)
save("data-mod-4.jld", "t", t, "x", x)
# 5 Hz
@time t, x = trapsolver(1,16,0.01,5)
save("data-mod-5.jld", "t", t, "x", x)
# 10 Hz
@time t, x = trapsolver(1,16,0.01,10)
save("data-mod-6.jld", "t", t, "x", x)
# 20 Hz
@time t, x = trapsolver(1,16,0.01,20)
save("data-mod-7.jld", "t", t, "x", x)
# 40 Hz
@time t, x = trapsolver(1,16,0.01,40)
save("data-mod-8.jld", "t", t, "x", x)


save("data-mod.jld", "t", t, "x", x)