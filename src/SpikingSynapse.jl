immutable SpikingSynapseParameters
  τpre::Float
  τpost::Float
  Apre::Float
  Apost::Float
end

export SpikingSynapse
type SpikingSynapse
  params::SpikingSynapseParameters
  rowptr::Array{Int,1}
  colptr::Array{Int,1}
  I::Array{Int,1}
  J::Array{Int,1}
  W::Array{Float,1}
  tpre::Array{Float,1}
  tpost::Array{Float,1}
  apre::Array{Float,1}
  apost::Array{Float,1}
  fireI::Array{Bool,1}
  fireJ::Array{Bool,1}
  g::Array{Float,1}
  records::Dict
end
function SpikingSynapse(pre, post, sym;
  σ=0.0, p=0.0, τpre=20, τpost=20, Apre=0.01, Apost=-0.01*1.05)
    params = SpikingSynapseParameters(τpre,τpost,Apre,Apost)
    w = σ*sprand(post.N,pre.N,p)
    rowptr, colptr, I, J, W = dsparse(w)
    tpre = zeros(W)
    tpost = zeros(W)
    apre = zeros(W)
    apost = zeros(W)
    fireI = post.fire
    fireJ = pre.fire
    g = getfield(post,sym)
    SpikingSynapse(params, rowptr, colptr, I, J, W, tpre, tpost, apre, apost, fireI, fireJ, g, Dict())
end

export forward!
@replace function forward!(c::SpikingSynapse, params::SpikingSynapseParameters, t)
  @inbounds for j in 1:length(colptr)-1
    if fireJ[j]
      for s in colptr[j]: (colptr[j+1] - 1)
        g[I[s]] += W[s]
      end
    end
  end
end


export integrate!
@replace function integrate!(c::SpikingSynapse, params::SpikingSynapseParameters, t)
end

export plasticity!
@replace function plasticity!(c::SpikingSynapse, params::SpikingSynapseParameters, t)
  @inbounds for j in 1:length(colptr)-1
    if fireJ[j]
      for s in colptr[j] : (colptr[j+1] - 1)
        apre[s] *= exp(- dt * (t - tpre[s]) / τpre)
        apost[s] *= exp(- dt * (t - tpost[s]) / τpost)
        apre[s] += Apre
        tpre[s] = t
        W[s]  += apost[s]
      end
    end
  end
  @inbounds for i in 1:length(rowptr)-1
    if fireI[i]
      for s in rowptr[i]: (rowptr[i+1] - 1)
        apre[s] *= exp(- dt * (t - tpre[s]) / τpre)
        apost[s] *= exp(- dt * (t - tpost[s]) / τpost)
        apost[s] += Apost
        tpost[s] = t
        W[s]  += apre[s]
      end
    end
  end
end


# immutable SpikingSynapseParameters
#   τpre::Float
#   τpost::Float
#   Apre::Float
#   Apost::Float
# end
#
# export SpikingSynapse
# type SpikingSynapse
#   params::SpikingSynapseParameters
#   rowptr::Array{Int,1}
#   colptr::Array{Int,1}
#   I::Array{Int,1}
#   J::Array{Int,1}
#   W::Array{Float,1}
#   apre::Array{Float,1}
#   apost::Array{Float,1}
#   fireI::Array{Bool,1}
#   fireJ::Array{Bool,1}
#   g::Array{Float,1}
#   records::Dict
# end
# function SpikingSynapse(pre, post, sym;
#   σ=0.0, p=0.0, τpre=20, τpost=20, Apre=0.01, Apost=-0.01*1.05)
#     params = SpikingSynapseParameters(τpre,τpost,Apre,Apost)
#     w = σ*sprand(post.N,pre.N,p)
#     rowptr, colptr, I, J, W = dsparse(w)
#     apre = zeros(W)
#     apost = zeros(W)
#     fireI = post.fire
#     fireJ = pre.fire
#     g = getfield(post,sym)
#     SpikingSynapse(params, rowptr, colptr, I, J, W, apre, apost, fireI, fireJ, g, Dict())
# end
#
# export forward!
# @replace function forward!(c::SpikingSynapse, params::SpikingSynapseParameters)
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
# @replace function integrate!(c::SpikingSynapse, params::SpikingSynapseParameters)
#   @inbounds for s in eachindex(apre)
#     apre[s] *= 1 - dt / τpre
#     apost[s] *= 1 - dt / τpost
#   end
# end
#
# export plasticity!
# @replace function plasticity!(c::SpikingSynapse, params::SpikingSynapseParameters)
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
