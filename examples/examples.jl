### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ e3cdb650-2d89-11f0-112e-c9bdfd106538
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 430ea2f1-2dd3-48c5-93ea-775894289b70
using PlutoLinks: @revise

# ╔═╡ 4282f8c6-568e-492d-83f3-b9b4e1c80392
using Plots

# ╔═╡ 43cdb6a1-0883-454b-bd9b-b7f9b3d22513
@revise using ParametricAdaptiveSampling

# ╔═╡ d358d19f-8e68-4c33-9936-cb5d16a1b74f
begin
	function plot_with_markers(plt, ps; plot_markers=true, kwargs...)
		plot!(plt, ps; label="f", kwargs...)
		if plot_markers
			@info "$(length(ps)) points"
			scatter!(plt, ps, markersize=1, makershape=:circle; label="samples")
		end
	end
	
	plot_with_markers(ps; kwargs...) = plot_with_markers(plot(), ps; kwargs...)
end

# ╔═╡ 18569ebf-448f-40e6-9e7f-53b304b7306a
"""
Plot first the curve (so you can see stuff), then the curve with markers, then a histogram of t values.
"""
function plot_parametric_sampling(ps_out; kwargs...)
	plts = plot(layout=(3,1), size=(600,1200))
	plot_with_markers(plts[1], ps_out[2]; plot_markers=false, kwargs...)
	plot_with_markers(plts[2], ps_out[2]; kwargs...)
	histogram!(plts[3], ps_out[1]; bins=:auto, label="t", ylabel="freq")
	plts
end

# ╔═╡ e6cb2d5a-04e7-4778-949f-fd345d42640f
md"""
# Some parametric examples
"""

# ╔═╡ e157d10d-a50f-410a-a0f9-5d17d483d61b
let f = t -> t .* (cos(t), sin(t))
	plot_parametric_sampling(sample_adaptive_parametric(f, 0, 4*pi))
end

# ╔═╡ 263da49d-9cfe-4405-8718-06efe39421a8
let f = t -> 1/t .* (cos(t), sin(t))
	plot_parametric_sampling(sample_adaptive_parametric(f, 1, 10*pi))
end

# ╔═╡ 449c99ac-9980-487e-b103-cfc3b7905346
let f = t -> (1+t*sin(5*t)) .* (cos(t), sin(t))
	plot_parametric_sampling(sample_adaptive_parametric(f, 1, 10*pi))
end

# ╔═╡ a276add9-c720-446a-b2ee-78d91cece448
md"""
# Some non-parametric examples
"""

# ╔═╡ 9422be81-223d-4bf7-9bbb-5dff926a9e16
plot_with_markers(sample_adaptive(x -> x^2, -10, 10))

# ╔═╡ a48d67a6-b4e2-4c38-aec9-d84a967e74e8
plot_with_markers(sample_adaptive(sin, 0, 2*pi))

# ╔═╡ fc4c7d06-09ed-40d8-a60e-541fcdbd2be4
# Graceful degredation of precision.
plot_with_markers(sample_adaptive(sin, 0, 2*pi; max_points=20))

# ╔═╡ 41eacfc6-8215-495a-8a9d-f56261012318
md"""
# Stress Tests
"""

# ╔═╡ 333f3810-bee8-4a70-8810-0ae9cf5411ee
let f = t -> (1/t, sin(t))
	plot_parametric_sampling(sample_adaptive_parametric(f, 1, 200))
end
# This breaches max_points, which it has every right to do.
# Note the pretty even distribution of t values. This is because two effects cancel each other out:
# - for high t, we don't have crazy variability in y direction but there's actually a lot of these in x direction in value space.
# - for low t, the other way round.

# ╔═╡ dd1ef647-80df-4b07-a9ea-e6cd3696782c
let f = t -> (1/t, 1/t*sin(t))
	plot_parametric_sampling(sample_adaptive_parametric(f, 1, 200))
end

# ╔═╡ 4ecb1a8a-4ce1-4ec3-9bc1-dbe063e60c13
let
	# We use our spiral from above but now the input is a pretty tight gaussian.
	# Note that the path goes down the spiral then back again, but this is not visible.
	f(t) = 1/t .* (cos(t), sin(t))
	g(t) = 10pi * exp(-t^2 * 100) + 1
	h(t) = f(g(t))
	plot_parametric_sampling(sample_adaptive_parametric(h, -5, 5))
end

# ╔═╡ 7093501c-b1ea-4e7d-a5c6-c128c8f3dd2c
# Non-differentiable points are generally fine.
# (NB we choose a weird xmax to avoid alignment with the initial samples grid)
plot_with_markers(sample_adaptive(abs, -1, 1.325))

# ╔═╡ 51dab9de-7b47-4943-9007-f8c661c1e35c
# Currently no implemented feature: Poles don't work and are not detected.
plot_with_markers(sample_adaptive(tan, 0, pi); ylim=(-20,20))

# ╔═╡ adbf7481-f15f-494f-b3f5-625a8ce69f18
md"""
# 3D Plots

This is supported just fine as well, as is any other number of dimensions. Requires the parametric version, though. (non-parametric 3d single-dimensional line plots are weird and not supported)

Note: the adaptive refinement does not consider perspective (of course). This tends to over-sample a bit.
"""

# ╔═╡ 2a28c543-c115-4f88-93d4-cb89e7eec67c
"""
[(a1, b1, c1), (a2, b2, c2)] -> ([a1, a2], [b1, b2], [c1, c2])
"""
function split_coords(vs)
	n = length(vs[1])
	Tuple([[v[i] for v in vs] for i in 1:n])
end

# ╔═╡ 67c4f572-4d93-489c-86f9-8a596080ddf5
let
	res = sample_adaptive_parametric(t -> (t^0.2, 1/t*cos(t), 1/t*sin(t)), 2, 50pi, tol=1e-3)[2]
	@info "$(length(res)) points"
	plot_with_markers(res |> split_coords)
end

# ╔═╡ Cell order:
# ╠═e3cdb650-2d89-11f0-112e-c9bdfd106538
# ╠═430ea2f1-2dd3-48c5-93ea-775894289b70
# ╠═4282f8c6-568e-492d-83f3-b9b4e1c80392
# ╠═43cdb6a1-0883-454b-bd9b-b7f9b3d22513
# ╠═d358d19f-8e68-4c33-9936-cb5d16a1b74f
# ╠═18569ebf-448f-40e6-9e7f-53b304b7306a
# ╟─e6cb2d5a-04e7-4778-949f-fd345d42640f
# ╠═e157d10d-a50f-410a-a0f9-5d17d483d61b
# ╠═263da49d-9cfe-4405-8718-06efe39421a8
# ╠═449c99ac-9980-487e-b103-cfc3b7905346
# ╟─a276add9-c720-446a-b2ee-78d91cece448
# ╠═9422be81-223d-4bf7-9bbb-5dff926a9e16
# ╠═a48d67a6-b4e2-4c38-aec9-d84a967e74e8
# ╠═fc4c7d06-09ed-40d8-a60e-541fcdbd2be4
# ╟─41eacfc6-8215-495a-8a9d-f56261012318
# ╠═333f3810-bee8-4a70-8810-0ae9cf5411ee
# ╠═dd1ef647-80df-4b07-a9ea-e6cd3696782c
# ╠═4ecb1a8a-4ce1-4ec3-9bc1-dbe063e60c13
# ╠═7093501c-b1ea-4e7d-a5c6-c128c8f3dd2c
# ╠═51dab9de-7b47-4943-9007-f8c661c1e35c
# ╠═adbf7481-f15f-494f-b3f5-625a8ce69f18
# ╠═2a28c543-c115-4f88-93d4-cb89e7eec67c
# ╠═67c4f572-4d93-489c-86f9-8a596080ddf5
