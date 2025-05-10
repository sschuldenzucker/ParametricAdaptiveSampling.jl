module ParametricAdaptiveSampling

export adaptive_sample_parametric, adaptive_sample

using LinearAlgebra: norm
using DataStructures
using IterTools: partition  # partition in lockstep
using Test

"Interval of t values, with function values. Ensures t1 < t2 but we don't check this on construction. V must be a tuple type."
struct Segment{V}
    t1::Float64
    v1::V
    t_mid::Float64
    v_mid::V
    t2::Float64
    v2::V
end

function Segment(f, t1, v1, t2, v2)
    t_mid = 0.5 * (t1 + t2)
    v_mid = f(t_mid)
    Segment(t1, v1, t_mid, v_mid, t2, v2)
end

function Segment(f, t1, t2)
    Segment(f, t1, f(t1), t2, f(t2))
end

mutable struct Range
    # Implicit: All lengths are equal
    mins::Vector{Float64}
    maxs::Vector{Float64}
    # Cached.
    widths::Vector{Float64}
end

Range(mins, maxs) = Range(mins, maxs, [maxs[i] - mins[i] for i in eachindex(mins)])

function Range(vs)
    v0 = first(vs)
    mins = [minimum(v[i] for v in vs) for i in eachindex(v0)]
    maxs = [maximum(v[i] for v in vs) for i in eachindex(v0)]
    Range(mins, maxs)
end

function empty_range(n::Integer)
    Range(fill(0.0, n), fill(0.0, n), fill(0.0, n))
end

"All lengths must be the same and v must use linear indexing"
function push_range!(range::Range, v)
    for i = 1:length(v)
        range.mins[i] = min(range.mins[i], v[i])
        range.maxs[i] = max(range.maxs[i], v[i])
        range.widths[i] = range.maxs[i] - range.mins[i]
    end
end

"Check error value against tolerance and push if above."
function maybe_push_errfun!(queue, errfun, tol, x, args...)
    err = errfun(x, args...)
    if err > tol
        push!(queue, (-err, x))
    end
end

"Error relative to the range. Component-wise."
function err_relative_range(s::Segment, range::Range)
    error_components = @. abs(s.v_mid - 0.5 * (s.v1 + s.v2))
    maximum(zip(error_components, range.widths)) do (err, w)
        # TODO handle 0. Is this good?
        err / max(w, eps())
    end
end

@test err_relative_range(
    Segment(t -> (t, t^2), 0.0, 10.0),
    Range([0.0, 0.0], [1000.0, 1000.0]),
) == 0.025

"""
    adaptive_sample_parametric(f, tmin, tmax)
    adaptive_sample_parametric(f, ts_init)

### Arguments
- `f`: The parameteric function. Given a Float64 in `[tmin, tmax]` must return an n-element tuple, where n must be constant across t values.
- `tmin::Real`: The lower bound of the parameter domain.
- `tmax::Real`: The upper bound of the parameter domain.
- `ts_init::Collection{Real}`: If this form is used, the initial t values are given explicitly. Otherwise (with `tmin` and `tmax`) we do initial sampling according to `min_points`. Must have size at least 2.

`f` must have well-defined values within `[tmin, tmax]`, including the endpoints.

### Keyword Arguments
- `errfun`: Error function `(Segment, ra::Range) -> Float64` where `ra` is (an estimate of) the value range. Default = l-infinity distance component-wise relative to range; this is a good default for plotting.
- `tol`: Tolerance to `errfun` for refinement. For plotting, this could be 1 / (max vertical/horizontal resolution). Default = 1e-3.
- `min_points::Integer`: Minimum number of t values to sample linearly initially. Pass 2 for no presampling beyond tmin and tmax.
- `max_points::Real`: Maximum number of points to sample. Pass `Inf` for no limit.

### Notes
- If `min_points` is low, `range` in `errfun` may be catastrophically small, at least for initial refinement. You probably don't want that.

### Features not yet implemented (SOMEDAY)
- Infinite or open t ranges where we would adaptively sample larger/smaller (towards 0) t values. This may be an instance of a more general hook into refinement.
- Optional penalty for recursion depth (or something like this) to avoid excessive concentration at polar points if `max_points` gets exhausted.
- Some way of detecting and avoiding polar points. Basically a way to say "it's not useful to plot this, better to cut it out".
- Option to drop points that are ultimately not needed if the range expands during refinement.
- Optionally, is there some smartness we can do with the derivative of f (calculated using ForwardDiff)?
"""
function adaptive_sample_parametric(
    f,
    tmin::Real,
    tmax::Real;
    tol = 1e-4,
    errfun = err_relative_range,
    min_points::Integer = 10,
    max_points = 4_000,
)
    # Initial uniform sampling
    @assert min_points >= 2
    ts_init = range(tmin, tmax, length = min_points)
    adaptive_sample_parametric(
        f,
        ts_init;
        tol = tol,
        errfun = errfun,
        max_points = max_points,
    )
end

function adaptive_sample_parametric(
    f,
    ts_init;
    tol = 1e-4,
    errfun = err_relative_range,
    max_points = 4_000,
)
    @assert length(ts_init) >= 2
    vs_init = f.(ts_init)

    # We store initial points here and append to it as we go (out of order, sorted later)
    # SOMEDAY smarter accounting to avoid the final sort?
    ps = zip(ts_init, vs_init) |> collect

    # Initialize range. We'll update as we go.
    ra = Range(vs_init)

    # Queue of pending splits, in error order, so we can stick to `max_points`. We sort again later.
    # Initialize with pairs of successive points.
    queue = BinaryHeap{Tuple{Float64,Segment}}(Base.By(first), [])
    for ((t1, v1), (t2, v2)) in partition(ps, 2, 1)
        s = Segment(f, t1, v1, t2, v2)
        maybe_push_errfun!(queue, errfun, tol, s, ra)
    end

    # Refinement.
    while !isempty(queue) && length(ps) < max_points
        (_, s) = pop!(queue)
        push!(ps, (s.t_mid, s.v_mid))
        s1 = Segment(f, s.t1, s.v1, s.t_mid, s.v_mid)
        s2 = Segment(f, s.t_mid, s.v_mid, s.t2, s.v2)
        push_range!(ra, s1.v_mid)
        push_range!(ra, s2.v_mid)

        maybe_push_errfun!(queue, errfun, tol, s1, ra)
        maybe_push_errfun!(queue, errfun, tol, s2, ra)
    end

    if !isempty(queue)
        max_err = first(queue)[1]
        @warn "max_points=$max_points reached with errors above tolerarnce. Max error = $(-max_err) > $tol = tolerance"
    end

    # SOMEDAY can we eliminate sorting by smarter accounting?
    sort!(ps; by = first)
    return map(first, ps), [p[2] for p in ps]
end

"""
    adaptive_sample(f, xmin, xmax)
    adaptive_sample(f, xs)

Non-parametric variant of `adaptive_sample_parametric()` for a function f: Number -> Number. See the documentation of that function.

Note: This function is often not needed. In many cases, for plotting, it will be faster to just spam a lot of x values instead. Cases where you _may_ want this function could be:

1. f is slow to evaluate.
2. You want to limit the number of data points for some later computation step that is slow.
3. You need very high precision in your output, maybe for a non-plotting use cases.

In cases 2. and 3., you may also want to use custom values and/or a custom error function in the keyword arguments.
"""
function adaptive_sample(
    f,
    xmin::Real,
    xmax::Real;
    tol = 1e-4,
    errfun = err_relative_range,
    min_points::Integer = 10,
    max_points = 4_000,
)
    # Initial uniform sampling
    @assert min_points >= 2
    xs_init = range(xmin, xmax, length = min_points)
    adaptive_sample(f, xs_init; tol = tol, errfun = errfun, max_points = max_points)
end

function adaptive_sample(f, xs; tol = 1e-4, errfun = err_relative_range, max_points = 4_000)
    f_parametric(x) = (x, f(x))
    adaptive_sample_parametric(
        f_parametric,
        xs,
        tol = tol,
        errfun = errfun,
        max_points = max_points,
    )
end

end # module ParametricAdaptiveSampling
