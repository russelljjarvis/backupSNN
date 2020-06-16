@snn_kw struct SpikingSynapseParameter
    τpre::Float32 = 20ms
    τpost::Float32 = 20ms
    Wmax::Float32 = 0.01
    ΔApre::Float32 = 0.01 * Wmax
    ΔApost::Float32 = -ΔApre * τpre / τpost * 1.05
end

@snn_kw mutable struct SpikingSynapse
    param::SpikingSynapseParameter = SpikingSynapseParameter()
    rowptr::Vector{Int32} # row pointer of sparse W
    colptr::Vector{Int32} # column pointer of sparse W
    I::Vector{Int32}      # postsynaptic index of W
    J::Vector{Int32}      # presynaptic index of W
    index::Vector{Int32}  # index mapping: W[index[i]] = Wt[i], Wt = sparse(dense(W)')
    W::Vector{Float32}  # synaptic weight
    tpre::Vector{Float32} = zero(W) # presynaptic spiking time
    tpost::Vector{Float32} = zero(W) # postsynaptic spiking time
    Apre::Vector{Float32} = zero(W) # presynaptic trace
    Apost::Vector{Float32} = zero(W) # postsynaptic trace
    fireI::Vector{Bool} # postsynaptic firing
    fireJ::Vector{Bool} # presynaptic firing
    g::Vector{Float32} # postsynaptic conductance
    records::Dict = Dict()
end

function SpikingSynapse(pre, post, sym; σ = 0.0, p = 0.0)
    w = σ * sprand(post.N, pre.N, p)
    rowptr, colptr, I, J, index, W = dsparse(w)
    fireI, fireJ = post.fire, pre.fire
    g = getfield(post, sym)
    SpikingSynapse(;@symdict(rowptr, colptr, I, J, index, W, fireI, fireJ, g)...)
end

function forward!(c::SpikingSynapse, param::SpikingSynapseParameter)
    @unpack colptr, I, W, fireJ, g = c
    @inbounds for j in 1:(length(colptr) - 1)
        if fireJ[j]
            for s in colptr[j]:(colptr[j+1] - 1)
                g[I[s]] += W[s]
            end
        end
    end
end

function plasticity!(c::SpikingSynapse, param::SpikingSynapseParameter, dt::Float32, t::Float32)
    @unpack rowptr, colptr, I, J, index, W, tpre, tpost, Apre, Apost, fireI, fireJ, g = c
    @unpack τpre, τpost, Wmax, ΔApre, ΔApost = param
    @inbounds for j in 1:(length(colptr) - 1)
        if fireJ[j]
            for s in colptr[j]:(colptr[j+1] - 1)
                Apre[s] *= exp32(- (t - tpre[s]) / τpre)
                Apost[s] *= exp32(- (t - tpost[s]) / τpost)
                Apre[s] += ΔApre
                tpre[s] = t
                W[s] = clamp(W[s] + Apost[s], 0f0, Wmax)
            end
        end
    end
    @inbounds for i in 1:(length(rowptr) - 1)
        if fireI[i]
            for st in rowptr[i]:(rowptr[i+1] - 1)
                s = index[st]
                Apre[s] *= exp32(- (t - tpre[s]) / τpre)
                Apost[s] *= exp32(- (t - tpost[s]) / τpost)
                Apost[s] += ΔApost
                tpost[s] = t
                W[s] = clamp(W[s] + Apre[s], 0f0, Wmax)
            end
        end
    end
end
