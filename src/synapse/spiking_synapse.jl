@with_kw struct SpikingSynapseParameter
    τpre::SNNFloat = 20ms
    τpost::SNNFloat = 20ms
    Wmax::SNNFloat = 0.01
    ΔApre::SNNFloat = 0.01 * Wmax
    ΔApost::SNNFloat = -ΔApre * τpre / τpost * 1.05
end

@with_kw mutable struct SpikingSynapse
    param::SpikingSynapseParameter = SpikingSynapseParameter()
    rowptr::Vector{SNNInt} # row pointer of sparse W
    colptr::Vector{SNNInt} # column pointer of sparse W
    I::Vector{SNNInt}      # postsynaptic index of W
    J::Vector{SNNInt}      # presynaptic index of W
    index::Vector{SNNInt}  # index mapping: W[index[i]] = Wt[i], Wt = sparse(dense(W)')
    W::Vector{SNNFloat}  # synaptic weight
    tpre::Vector{SNNFloat} = zero(W) # presynaptic spiking time
    tpost::Vector{SNNFloat} = zero(W) # postsynaptic spiking time
    Apre::Vector{SNNFloat} = zero(W) # presynaptic trace
    Apost::Vector{SNNFloat} = zero(W) # postsynaptic trace
    fireI::Vector{Bool} # postsynaptic firing
    fireJ::Vector{Bool} # presynaptic firing
    g::Vector{SNNFloat} # postsynaptic conductance
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

function plasticity!(c::SpikingSynapse, param::SpikingSynapseParameter, dt::SNNFloat, t::SNNFloat)
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
