"""
    A macro equivalent to @batch per=thread from Polyester.jl.
"""
macro tbatch(ex)
    reserve, minbatch, per, threadlocal = 0, 1, :thread, (Symbol(""), :Any)
    Polyester.enclose(macroexpand(__module__, ex), reserve, minbatch, per, threadlocal, __module__)
end
