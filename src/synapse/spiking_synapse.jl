@withkw immutable SpikingSynapseParameter
    τpre::Float = 20ms
    τpost::Float = 20ms
    Wmax::Float = 0.01
    ΔApre::Float = 0.01 * Wmax
    ΔApost::Float = -ΔApre * τpre / τpost * 1.05
end

@withkw type SpikingSynapse
    param::SpikingSynapseParameter = SpikingSynapseParameter()
    rowptr::Vector{Int} # row pointer of sparse W
    colptr::Vector{Int} # column pointer of sparse W
    I::Vector{Int}      # postsynaptic index of W
    J::Vector{Int}      # presynaptic index of W
    index::Vector{Int}  # index mapping: W[index[i]] = Wt[i], Wt = sparse(dense(W)')
    W::Vector{Float}  # synaptic weight
    tpre::Vector{Float} = zeros(W) # presynaptic spiking time
    tpost::Vector{Float} = zeros(W) # postsynaptic spiking time
    Apre::Vector{Float} = zeros(W) # presynaptic trace
    Apost::Vector{Float} = zeros(W) # postsynaptic trace
    fireI::Vector{Bool} # postsynaptic firing
    fireJ::Vector{Bool} # presynaptic firing
    g::Vector{Float} # postsynaptic conductance
    records::Dict = Dict()
end

function SpikingSynapse(pre, post, sym; σ = 0.0, p = 0.0)
    w = σ * sprand(post.N, pre.N, p)
    rowptr, colptr, I, J, index, W = dsparse(w)
    fireI, fireJ = post.fire, pre.fire
    g = getfield(post, sym)
    SpikingSynapse(;@symdict(rowptr, colptr, I, J, index, W, fireI, fireJ, g)...)
end

@replace function forward!(c::SpikingSynapse, param::SpikingSynapseParameter)
    @inbounds for j in 1:(length(colptr) - 1)
        if fireJ[j]
            for s in colptr[j]:(colptr[j+1] - 1)
                g[I[s]] += W[s]
            end
        end
    end
end

@replace function plasticity!(c::SpikingSynapse, param::SpikingSynapseParameter, dt::Float, t::Float)
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
