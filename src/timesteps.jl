"""
    Some code to simplify the use of several time steps that progress 
    simultaneously in a loop.
"""
mutable struct TimeStepper{T}
    i::Int
    start::T
    step::T
end

function TimeStepper(start, step)
    TimeStepper(0, start, step)
end

function TimeStepper(step)
    TimeStepper(zero(step), step)
end

""" Checks if it is time to advance the counter.  In that case first calls
    function `f` as `f(i)` where `i` is the counter _before_ the increase.
"""
function atstep(f::Function, ts::TimeStepper, t)
    if t >= ts.step * ts.i
        f(ts.i)
        ts.i += 1
    end
end
 
