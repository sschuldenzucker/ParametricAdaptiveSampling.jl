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
    factor = 0.5
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

"Push to a queue and cache the error value. `errfun: x |-> Number`"
function push_errfun!(queue, errfun, x, args...)
    push!(queue, (-errfun(x, args...), x))
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
- `tmin::Float64`: The lower bound of the parameter domain.
- `tmax::Float64`: The upper bound of the parameter domain.

### Keyword Arguments
- `errfun`: Error function `(Segment, range::Tuple{Float64,...}) -> Float64` where `range` is (an estimate of) the value range. Default = relative to range; this is a good default for plotting.
- `tolerance`: Tolerance for refinement.
- `min_points::Int`: Minimum number of points to sample initially.
- `max_points::Int` or `Nothing`: Maximum number of points to sample. Should be `> min_points`. `nothing` means unlimited.

### Notes
- If `min_points` is low, `range` in `errfun` may be catastrophically imprecise, at least for initial refinement. You probably don't want that.

### Features not yet implemented (SOMEDAY)
- Allow to give the initial t values directly (to vary the presampling strategy)
- Allow to configure the refinement strategy in a more detailed way (see Notes).
- Infinite or open t ranges where we would adaptively sample larger/smaller (towards 0) t values.
- Drop points that are ultimately not needed if the range expands during refinement.
"""
function adaptive_sample_parametric(
    f::Function,
    tmin::Float64,
    tmax::Float64;
    tolerance = 1e-3,
    errfun = err_relative_range,
    min_points = 10,
    max_points = 4_000,
)
    # todo min_points, max_points option for not given. Can we use Inf for max and 2 for min?

    # Queue in error order, so we can stick to `max_points`. We sort again later.
    # SOMEDAY can we eliminate that by smarter accounting?
    queue = BinaryHeap{Tuple{Float64,Segment}}(Base.By(first), [])

    # Initial uniform sampling
    @assert min_points >= 2
    # SOMEDAY may want to sample a single value first to determine the value type (specifically, )
    ts_init = range(tmin, tmax, length = min_points)
    vs_init = f.(ts_init)
    # successive pairs of initial points
    pairs_init = partition(zip(ts_init, vs_init), 2, 1)
    segments_init = [Segment(f, t1, v1, t2, v2) for ((t1, v1), (t2, v2)) in pairs_init]

    # We use the value at tmax to presample and to understand the dimension.
    # (TODO should we?)

    # Fill this from initial segments. We'll update as we go.
    ra = Range(vs_init)
    for s in segments_init
        push_errfun!(queue, errfun, s, ra)
    end
    @debug queue

    # Refinement and commit.
    committed = []
    # TODO max_points not working well, getting stuck in infinite loop. Should consider length(queue) as well.
    # But if we get too many points *this* way, it's kinda sad: alg is clearly stuck.
    # This could be b/c of (1) polar points and (2) values close to 0.
    # Maybe we should commit earlier?
    while !isempty(queue) && length(committed) < max_points
        (err, s) = pop!(queue)
        if err <= tolerance
            push!(committed, (s.t1, s.v1))
        else
            s1 = Segment(f, s.t1, s.v1, s.t_mid, s.v_mid)
            s2 = Segment(f, s.t_mid, s.v_mid, s.t2, s.v2)
            push_range!(ra, s1.v_mid)
            push_range!(ra, s2.v_mid)
            push_errfun!(queue, errfun, s1, ra)
            push_errfun!(queue, errfun, s2, ra)
        end
    end
    # Accounting not quite perfect here.
    push!(committed, (tmax, f(tmax)))

    sort!(committed; by = first)

    return map(first, committed), [p[2] for p in committed]
end

end # module ParametricAdaptiveSampling
