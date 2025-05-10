module ParametricAdaptiveSampling

export adaptive_sample_parametric

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

function Segment(f::Function, t1, v1, t2, v2)
    t_mid = 0.5 * (t1 + t2)
    v_mid = f(t_mid)
    # abs_error = norm(@. v_mid - 0.5 * (v1 + v2))
    # # SOMEDAY we should use some sane value instead of eps I think. Allow v_mid to be 0.
    # rel_error = abs_error / max(norm(v_mid), eps())
    Segment(t1, v1, t_mid, v_mid, t2, v2)
end

function Segment(f::Function, t1, t2)
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

# TODO review type stability, maybe make a function accepting v_tmax as a sample, and maybe use static arrays to capture the length.

"""
### Arguments
- `f::Function`: The parameteric function. Given a Float64 in `[tmin, tmax]` must return an n-element tuple, where n must be constant across t values.
- `tmin::Real`: The lower bound of the parameter domain.
- `tmax::Real`: The upper bound of the parameter domain.

`f` must have well-defined values within `[tmin, tmax]`, including the endpoints.

### Keyword Arguments
- `errfun`: Error function `(Segment, ra::Range) -> Float64` where `ra` is (an estimate of) the value range. Default = l-infinity distance component-wise relative to range; this is a good default for plotting.
- `tol`: Tolerance to `errfun` for refinement. For plotting, this could be 1 / (max vertical/horizontal resolution). Default = 1e-3.
- `min_points::Integer`: Minimum number of points to sample initially. Pass 2 for no presampling beyond tmin and tmax.
- `max_points::Real`: Maximum number of points to sample. Should be `>= min_points`. Pass `Inf` for no limit.

### Notes
- If `min_points` is low, `range` in `errfun` may be catastrophically small, at least for initial refinement. You probably don't want that.

### Features not yet implemented (SOMEDAY)
- Allow to give the initial t values directly, rather than an implicit linear range (to vary the presampling strategy)
- Infinite or open t ranges where we would adaptively sample larger/smaller (towards 0) t values. This may be an instance of a more general hook into refinement.
- Optional penalty for recursion depth (or something like this) to avoid excessive concentration at polar points if `max_points` gets exhausted.
- Option to drop points that are ultimately not needed if the range expands during refinement.
"""
function adaptive_sample_parametric(
    f::Function,
    tmin::Real,
    tmax::Real;
    tol = 1e-4,
    errfun = err_relative_range,
    min_points::Integer = 10,
    max_points = 4_000,
)
    # Initial uniform sampling
    @assert min_points >= 2
    # SOMEDAY may want to sample a single value first to determine the value type (specifically, )
    ts_init = range(tmin, tmax, length = min_points)
    vs_init = f.(ts_init)

    # We store initial points here and append to it as we go (out of order, sorted later)
    # SOMEDAY smarter accounting to avoid the final sort?
    ps = zip(ts_init, vs_init) |> collect

    # Initialize range. We'll update as we go.
    ra = Range(vs_init)

    # We use the value at tmax to presample and to understand the dimension.
    # (TODO should we?)

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
        @warn "max_points=$max_points reached with errors above tolerarnce. Max error = $max_err > $tol = tolerance"
    end

    # SOMEDAY can we eliminate sorting by smarter accounting?
    sort!(ps; by = first)
    return map(first, ps), [p[2] for p in ps]
end

end # module ParametricAdaptiveSampling
