module ParametricAdaptiveSampling

export sample_adaptive_parametric, sample_adaptive
public Segment, Range, widths, err_relative_range

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

"A range of vectors, from the 'lower-left' corner to the 'upper-right' one."
mutable struct Range
    # Implicit: All lengths are equal
    mins::Vector{Float64}
    maxs::Vector{Float64}
end

"Iterator for the widths of each dimension. An iterator of Float64."
widths(ra::Range) = (ma-mi for (ma, mi) in zip(ra.maxs, ra.mins))

function Range(vs)
    v0 = first(vs)
    mins = [minimum(v[i] for v in vs) for i in eachindex(v0)]
    maxs = [maximum(v[i] for v in vs) for i in eachindex(v0)]
    Range(mins, maxs)
end

function empty_range(n::Integer)
    Range(fill(0.0, n), fill(0.0, n))
end

"All lengths must be the same and v must use linear indexing"
function push_range!(range::Range, v)
    for i = 1:length(v)
        range.mins[i] = min(range.mins[i], v[i])
        range.maxs[i] = max(range.maxs[i], v[i])
    end
end

"Check error value against tolerance and push if above."
function maybe_push_errfun!(queue, errfun, tol, x, args...)
    err = errfun(x, args...)
    if err > tol
        push!(queue, (-err, x))
    end
end

"An `errfun` for the sampling functions. l-infinity error relative to the range, component-wise."
function err_relative_range(s::Segment, ra::Range)
    error_components = @. abs(s.v_mid - 0.5 * (s.v1 + s.v2))
    maximum(zip(error_components, widths(ra))) do (err, w)
        # TODO handle 0. Is this good? Maybe we should just raise a descriptive error.
        err / max(w, eps())
    end
end

@test err_relative_range(
    Segment(t -> (t, t^2), 0.0, 10.0),
    Range([0.0, 0.0], [1000.0, 1000.0]),
) == 0.025

"""
    sample_adaptive_parametric(f, tmin, tmax)
    sample_adaptive_parametric(f, ts_init)

### Arguments
- `f`: The parameteric function. Given a Float64 in `[tmin, tmax]` must return an n-element tuple or vector of reals, where n must be constant across t values. This should be at least continuous on the interval.
- `tmin::Real`: The lower bound of the parameter domain.
- `tmax::Real`: The upper bound of the parameter domain.
- `ts_init::Collection{Real}`: If this form is used, the initial t values are given explicitly. Otherwise (with `tmin` and `tmax`) we do initial sampling according to `min_points`. Must have size at least 2.

`f` must have well-defined values within `[tmin, tmax]`, including the endpoints.

### Keyword Arguments
- `errfun`: Error function `(Segment, ra::Range) -> Float64` where `ra` is (an estimate of) the value range. Default = l-infinity distance component-wise relative to range; this is a good default for plotting.
- `tol`: Tolerance to `errfun` for refinement. For plotting, this could be 1 / (max vertical/horizontal resolution). Default = 1e-3.
- `min_points::Integer`: Minimum number of t values to sample linearly initially. Pass 2 for no presampling beyond tmin and tmax.
- `max_points::Real`: Maximum number of points to sample. Pass `Inf` for no limit.

### Returns

Returns `(t values::Vector{Number}, value points::Vector{Tuple})`. For plotting applications, you'll want to discard the t values, i.e., take the second return value only.

### Notes

This works as follows: first, it performs an initial presampling step (first form) or accepts a list of initial samples (second form). Then, it considers the segment where the error of linear interpolation across the segment vs the midpoint value (midpoint in t space) is worst and splits it into two halves. Then repeat until all errors are below the tolerance or we run out of points. The estimate of the value range is updated as we add new points.

### Caveats

If `min_points` is low, `range` in `errfun` may be catastrophically small, at least for initial refinement. You probably don't want that. More generally, if the initial samples do not represent the value range reasonably well, in some cases, refinement may add excessive points, and may even reach `max_points` without any meaningful progress. If this happens, you probably want to adjust initial samples.

Recursion is is by splitting into halves, so if your function has details that are missed by this recursive grid, they won't show up. In this case, increasing `min_points` may help. 

### Features not yet implemented (SOMEDAY)

- Infinite or open t ranges where we would adaptively sample larger/smaller (towards infinity or the edges) t values. This may be an instance of a more general hook into refinement.
- Some way of detecting and avoiding polar points. Basically a way to say "it's not useful to plot this, better to cut it out". Not clear how we'd do this in a robust way without affecting some functions negatively.
- Optional penalty for recursion depth (or something like this) to avoid excessive concentration at polar points if `max_points` gets exhausted. This would have to be configured by the user with knowledge of the function, can be detrimental otherwise. Unclear if we want this.
- Option to split into more than two parts per recursion step. This could help when details are missed by the 2-split recursion (see Caveats).
- Point density target in value space. Could help in situations where the midway point is interpolated well, but there is additional variation at some other point.
- Option to drop points that are ultimately not needed if the range expands during refinement (may reduce number of points, might be useful for some later computation steps).
- Is there some smartness we can do with if the derivative of f is available?
- A way to understand that the `ra::Range` below may actually be the wrong thing, e.g., when we use
  `aspect_ratio=1` in the plot and the plot therefore creates additional space, or some xlim or
  ylim. Hard to do generically. Maybe make an option to pass the plot range explicitly.
"""
function sample_adaptive_parametric(
    f,
    tmin::Real,
    tmax::Real;
    tol = 5e-4,
    errfun = err_relative_range,
    min_points::Integer = 20,
    max_points = 4_000,
)
    # Initial uniform sampling
    @assert min_points >= 2
    ts_init = range(tmin, tmax, length = min_points)
    sample_adaptive_parametric(
        f,
        ts_init;
        tol = tol,
        errfun = errfun,
        max_points = max_points,
    )
end

function sample_adaptive_parametric(
    f,
    ts_init;
    tol = 5e-4,
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
        @warn "max_points=$max_points reached with errors above tolerance. Max error = $(-max_err) > $tol = tolerance"
    end

    # SOMEDAY can we eliminate sorting by smarter accounting?
    sort!(ps; by = first)
    # Convert to Tuple so plot() does the right thing. (kinda sad we have to do this)
    return map(first, ps), [Tuple(p[2]) for p in ps]
end

"""
    sample_adaptive(f, xmin, xmax)
    sample_adaptive(f, xs)

Non-parametric variant of `sample_adaptive_parametric()` for a function f: Number -> Number. See the documentation of that function.

### Returns

This returns a list of pairs (2-tuples), which are the generated points.

Note that the return format is different from `sample_adaptive_parametric()` because there is no separate space of "t values" here.

### When this function is useful

Note: This function is often not needed. In many cases, for plotting, it will be faster to just spam a lot of x values instead. Cases where you _may_ want this function could be:

1. f is slow to evaluate.
2. You want to limit the number of data points for some later computation step that is slow.
3. You need very high precision in your output and can't affort _that_ many points, maybe for a non-plotting use case.

In cases 2. and 3., you may also want to use custom values and/or a custom error function in the keyword arguments.
"""
function sample_adaptive(
    f,
    xmin::Real,
    xmax::Real;
    tol = 5e-4,
    errfun = err_relative_range,
    min_points::Integer = 10,
    max_points = 4_000,
)
    # Initial uniform sampling
    @assert min_points >= 2
    xs_init = range(xmin, xmax, length = min_points)
    sample_adaptive(f, xs_init; tol = tol, errfun = errfun, max_points = max_points)
end

function sample_adaptive(f, xs; tol = 5e-4, errfun = err_relative_range, max_points = 4_000)
    f_parametric(x) = (x, f(x))
    sample_adaptive_parametric(
        f_parametric,
        xs,
        tol = tol,
        errfun = errfun,
        max_points = max_points,
    )[2]
end

end # module ParametricAdaptiveSampling
