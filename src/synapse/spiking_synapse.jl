@withkw immutable SpikingSynapseParameter
  τpre::Float32 = 20
  τpost::Float32 = 20
  Apre::Float32 = 0.01
  Apost::Float32 = -0.01 * 1.05
end

@withkw type SpikingSynapse
  param::SpikingSynapseParameter = SpikingSynapseParameter()
  rowptr::Vector{Int}
  colptr::Vector{Int}
  I::Vector{Int}
  J::Vector{Int}
  W::Vector{Float32}
  tpre::Vector{Float32} = zeros(W)
  tpost::Vector{Float32} = zeros(W)
  apre::Vector{Float32} = zeros(W)
  apost::Vector{Float32} = zeros(W)
  fireI::Vector{Bool}
  fireJ::Vector{Bool}
  g::Vector{Float32}
  records::Dict = Dict()
end

function SpikingSynapse(pre, post, sym; σ = 0.0, p = 0.0)
  w = σ * sprand(post.N, pre.N, p)
  rowptr, colptr, I, J, W = dsparse(w)
  fireI, fireJ = post.fire, pre.fire
  g = getfield(post, sym)
  SpikingSynapse(;@symdict(rowptr, colptr, I, J, W, fireI, fireJ, g)...)
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

@replace function plasticity!(c::SpikingSynapse, param::SpikingSynapseParameter, dt::Float32, t::Float32)
  @inbounds for j in 1:(length(colptr) - 1)
    if fireJ[j]
      for s in colptr[j]:(colptr[j+1] - 1)
        apre[s] *= exp(- dt * (t - tpre[s]) / τpre)
        apost[s] *= exp(- dt * (t - tpost[s]) / τpost)
        apre[s] += Apre
        tpre[s] = t
        W[s]  += apost[s]
      end
    end
  end
  @inbounds for i in 1:(length(rowptr) - 1)
    if fireI[i]
      for s in rowptr[i]:(rowptr[i+1] - 1)
        apre[s] *= exp(- dt * (t - tpre[s]) / τpre)
        apost[s] *= exp(- dt * (t - tpost[s]) / τpost)
        apost[s] += Apost
        tpost[s] = t
        W[s]  += apre[s]
      end
    end
  end
end


# immutable SpikingSynapseParameter
#   τpre::Float32
#   τpost::Float32
#   Apre::Float32
#   Apost::Float32
# end
#
# export SpikingSynapse
# type SpikingSynapse
#   param::SpikingSynapseParameter
#   rowptr::Vector{Int}
#   colptr::Vector{Int}
#   I::Vector{Int}
#   J::Vector{Int}
#   W::Vector{Float32}
#   apre::Vector{Float32}
#   apost::Vector{Float32}
#   fireI::Vector{Bool}
#   fireJ::Vector{Bool}
#   g::Vector{Float32}
#   records::Dict
# end
# function SpikingSynapse(pre, post, sym;
#   σ=0.0, p=0.0, τpre=20, τpost=20, Apre=0.01, Apost=-0.01*1.05)
#     param = SpikingSynapseParameter(τpre,τpost,Apre,Apost)
#     w = σ*sprand(post.N,pre.N,p)
#     rowptr, colptr, I, J, W = dsparse(w)
#     apre = zeros(W)
#     apost = zeros(W)
#     fireI = post.fire
#     fireJ = pre.fire
#     g = getfield(post,sym)
#     SpikingSynapse(param, rowptr, colptr, I, J, W, apre, apost, fireI, fireJ, g, Dict())
# end
#
# export forward!
# @replace function forward!(c::SpikingSynapse, param::SpikingSynapseParameter)
#   @inbounds for j in 1:length(colptr)-1
#     if fireJ[j]
#       for s in colptr[j]:colptr[j+1]-1
#         g[I[s]] += W[s]
#       end
#     end
#   end
# end
#
#
# export integrate!
# @replace function integrate!(c::SpikingSynapse, param::SpikingSynapseParameter)
#   @inbounds for s in eachindex(apre)
#     apre[s] *= 1 - dt / τpre
#     apost[s] *= 1 - dt / τpost
#   end
# end
#
# export plasticity!
# @replace function plasticity!(c::SpikingSynapse, param::SpikingSynapseParameter)
#   @inbounds for j in 1:length(colptr)-1
#     if fireJ[j]
#       for s in colptr[j]:colptr[j+1]-1
#         apre[s] += Apre
#         W[s]  += apost[s]
#       end
#     end
#   end
#   @inbounds for i in 1:length(rowptr)-1
#     if fireI[i]
#       for s in rowptr[i]:rowptr[i+1]-1
#         apost[s] += Apost
#         W[s]  += apre[s]
#       end
#     end
#   end
# end
