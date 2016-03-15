__precompile__()
module SpikingNetworks
include("utility.jl")
include("codegen.jl")
include("group.jl")
include("unit.jl")
include("parse.jl")
include("init.jl")
include("stateMonitor.jl")
include("spike.jl")
include("reset.jl")
include("synapse.jl")
include("plasticity.jl")
end # module