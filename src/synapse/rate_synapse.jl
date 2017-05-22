@withkw immutable RateSynapseParameter
  lr::Float32 = 1e-3
end

@withkw type RateSynapse
  param::RateSynapseParameter = RateSynapseParameter()
  rowptr::Vector{Int}
  colptr::Vector{Int}
  I::Vector{Int}
  J::Vector{Int}
  W::Vector{Float32}
  rI::Vector{Float32}
  rJ::Vector{Float32}
  g::Vector{Float32}
  records::Dict = Dict()
end

function RateSynapse(pre, post; σ = 0.0, p = 0.0)
  w = σ / √(p * pre.N) * sprandn(post.N, pre.N, p)
  rowptr, colptr, I, J, W = dsparse(w)
  rI, rJ, g = post.r, pre.r, post.g
  RateSynapse(;@symdict(rowptr, colptr, I, J, W, rI, rJ, g)...)
end

@replace function forward!(c::RateSynapse, param::RateSynapseParameter)
    fill!(g, zero(eltype(g)))
    @inbounds for j in 1:(length(colptr) - 1)
      rJj = rJ[j]
       for s = colptr[j]:(colptr[j+1] - 1)
          g[I[s]] += W[s] * rJj
      end
    end
end

@replace function plasticity!(c::RateSynapse, param::RateSynapseParameter, dt::Float32, t::Float32)
  @inbounds for j in 1:(length(colptr) - 1)
    rowind = colptr[j]:(colptr[j+1] - 1)
    rIW = zero(Float32)
    for s in rowind
      rIW += rI[I[s]] * W[s]
    end
    Δ = lr*(rJ[j] - rIW)
    for s in rowind
        W[s] += rI[I[s]] * Δ
    end
  end
end
