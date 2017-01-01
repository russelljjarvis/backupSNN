immutable RateSynapseParameters
  lr::Float
end

export RateSynapse
type RateSynapse
  params::RateSynapseParameters
  rowptr::Array{Int,1}
  colptr::Array{Int,1}
  I::Array{Int,1}
  J::Array{Int,1}
  W::Array{Float,1}
  rI::Array{Float,1}
  rJ::Array{Float,1}
  g::Array{Float,1}
  records::Dict
end
function RateSynapse(pre, post; lr=1e-3, σ=0.0, p=0.0)
    params = RateSynapseParameters(lr)
    w = σ/√(p*pre.N)*sprandn(post.N, pre.N, p)
    rowptr, colptr, I, J, W = dsparse(w)
    rI = post.r
    rJ = pre.r
    g = post.g
    records = Dict()
    RateSynapse(params, rowptr, colptr, I, J, W, rI, rJ, g, records)
end

export forward!
@replace function forward!(c::RateSynapse, params::RateSynapseParameters, t)
    fill!(g, zero(eltype(g)))
    @inbounds for j in 1:length(colptr)-1
      rJj = rJ[j]
       for s = colptr[j] : (colptr[j+1]-1)
          g[I[s]] += W[s]*rJj
      end
    end
end

export integrate!
function integrate!(c::RateSynapse, params::RateSynapseParameters, t)
end

export plasticity!
@replace function plasticity!(c::RateSynapse, params::RateSynapseParameters, t)
  @inbounds for j in 1:length(colptr)-1
    rowind = colptr[j] : (colptr[j+1]-1)
    rIW = zero(Float)
    for s in rowind
      rIW += rI[I[s]] * W[s]
    end
    Δ = lr*(rJ[j] - rIW)
    for s in rowind
        W[s] += rI[I[s]] * Δ
    end
  end
end
