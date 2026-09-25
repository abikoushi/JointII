module empest

using LinearAlgebra
using Distributions
using Random
using LogExpFunctions
using SpecialFunctions
using Statistics
using SparseArrays


function prob_sample(x)
  out = [x[i] > 0.0 ? rand(Random.default_rng(), Gamma(x[i],1)) : 0.0 for i in eachindex(x)]
  out ./= sum(out)
  return out
end

function freq_sample(x)
  hs = [x[i] > 0.0 ? rand(Random.default_rng(), Gamma(x[i],1)) : 0.0 for i in eachindex(x)]
  return cumsum(hs)
end

#prob2ccdf(x) = 1 .- cumsum(x)
prob2ccdf(x) = reverse(cumsum(reverse(x)))

colmarginal(x) = cumsum(vec(sum(x,dims=2)))

function asample(le_rank,re_rank,ls_rank,rs_rank,h)
  m = size(h,1)
  a = spzeros(m,m)
  for i in (le_rank):(re_rank)
    start = max(ls_rank,i+1)
    for j in start:rs_rank
        a[i,j] += h[j-i]
    end
  end
  rind, cind, val = findnz(a)
  val ./= sum(val)
  occ = rand(Random.default_rng(), Categorical(val))
  out = zero(a)
  out[rind[occ], cind[occ]] = 1.0
  return out
end

function Asample!(a, le_rank,re_rank,ls_rank,rs_rank,h)
  m = size(h,1)
  A = zero(a)
  for i in eachindex(rs_rank)
    A .+= asample(le_rank[i], re_rank[i], ls_rank[i], rs_rank[i], h)
  end
  copy!(a, A)
end


function summingup!(α, A, h)
  lp = 0.
  m = size(A,1)
  alpha = zero(α)

  for i in 1:(m-1)
    for j in (i+1):m
      alpha[j-i] += A[i,j]
      if h[j-i]>0
      lp += xlogy(A[i,j],h[j-i])
      end
    end
  end
  copy!(α,alpha)
  return lp
end

function gibbssampler(LE, RE, LS, RS, iter)
  ti = sort(unique([LE;RE;LS;RS]))
  ti = ti[isfinite.(ti)]
  #n = size(LS,1)
  le_rank = indexin(LE, ti)
  re_rank = indexin(RE, ti)
  ls_rank = indexin(LS, ti)
  rs_rank = indexin(RS, ti)
  N = length(LE)
  m = size(ti,1)
  ccdf = zeros(iter, m)
  alpha = ones(m)
  A = spzeros(m, m)
  logprob = zeros(iter)
  for i in 1:iter
    h = prob_sample(alpha)
    Asample!(A, le_rank, re_rank, ls_rank, rs_rank, h)
    logprob[i] = summingup!(alpha, A, h)
    ccdf[i,:] = prob2ccdf(h)
  end
  return (value = ti, ccdf = ccdf, logprob = logprob)
end

end