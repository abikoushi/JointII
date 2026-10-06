module contool

using LinearAlgebra
using Distributions
using Random
using LogExpFunctions
using SpecialFunctions
using Statistics
using SparseArrays
using StatsBase

function expmeanlog(x)
  sumx = sum(x)
  [x[i]>0 ? exp(digamma(x[i]) - digamma(sumx)) : 0 for i in eachindex(x)]
end

#singly interval censored
function makeIC(rng::AbstractRNG, dt::UnivariateDistribution, N::Int)
  at = 0.
  LE = zeros(N)
  RE = zeros(N)
  S = zeros(N)
  for i in 1:N
      ue = rand(rng)
      at -= log(rand(rng))
      tau = -log(rand(rng))
      y = rand(rng, dt)
      S[i] = at + y
      LE[i] = max(0, at - tau * (1 - ue))
      RE[i] = min(S[i], at + tau * ue)
  end
  return LE, RE, S
end


function makeDIC(rng::AbstractRNG, dt::UnivariateDistribution, tau, WE, WS, shape=1, maxit=5000)
  at = 0.
  LE = zeros(maxit)
  RE = zeros(maxit)
  LS = zeros(maxit)
  RS = zeros(maxit)
  z = rand(rng)
  at = (-log(z))^(1/shape)
  for i in 1:maxit
      ue = rand(rng)
      us = rand(rng)
      y = rand(rng, dt)
      S = at + y
      LE[i] = max(0, at - WE * (1 - ue))
      LS[i] = max(LE[i], S - WS * (1 - us))
      RS[i] = S + WS * us
      RE[i] = min(RS[i], at + WE * ue)
      z = rand(rng) #next infection
      at =  ( (-log(z)) + at^shape)^(1/shape)
      if at >= tau
        LE=LE[1:i]
        LS=LS[1:i]
        RE=RE[1:i]
        RS=RS[1:i]
        break
      end
  end
  return LE, RE, LS, RS
end


function aup(le_rank,re_rank,ls_rank,rs_rank,h)
  m = size(h,1)
  a = spzeros(m,m)
  for i in (le_rank):(re_rank)
    start = max(ls_rank,i+1)
    for j in start:rs_rank
        a[i,j] += h[j-i]
        #print(i,", ",j,"\n")
    end
  end
  a ./= sum(a)
  return a
end

function Aup!(a, le_rank,re_rank,ls_rank,rs_rank,h)
  m = size(h,1)
  A = zero(a)
  for i in eachindex(rs_rank)
    A .+= aup(le_rank[i],re_rank[i],ls_rank[i],rs_rank[i],h)
  end
  copy!(a, A)
end




function summingup!(α, A, h)
  lp = 0.
  m = size(A,1)
  alpha = zero(α)
  for i in 1:(m-1)
    for j in (i+1):m
    #print(i,", ", j ,"\n")
      alpha[j-i] += A[i,j]
      if h[j-i]>0
      lp += xlogy(A[i,j],h[j-i])
      end
    end
  end
  copy!(α,alpha)
  return lp
end

function vb(LE, RE, LS, RS, maxit=1000, reltol=1e-8)
  ti = sort(unique([LE;RE;LS;RS]))
  ti = ti[isfinite.(ti)]
  #n = size(LS,1)
  le_rank = indexin(LE, ti)
  re_rank = indexin(RE, ti)
  ls_rank = indexin(LS, ti)
  rs_rank = indexin(RS, ti)
  m = size(ti,1)
  alpha = ones(m)
  A = spzeros(m,m)
  logprob = zeros(maxit)
  for i in 2:maxit
    h = expmeanlog(alpha)
    Aup!(A,le_rank, re_rank, ls_rank, rs_rank, h)
    logprob[i] = summingup!(alpha, A, h)
    if abs((logprob[i]-logprob[i-1])/logprob[i-1])<reltol
      logprob = logprob[2:i]
      break
    end
  end
  return (value=ti, alpha=alpha, A=A, logprob=logprob)
end

function EM(LE, RE, LS, RS, maxit=1000, reltol=1e-8)
  ti = sort(unique([LE;RE;LS;RS]))
  ti = ti[isfinite.(ti)]
  #n = size(LS,1)
  le_rank = indexin(LE, ti)
  re_rank = indexin(RE, ti)
  ls_rank = indexin(LS, ti)
  rs_rank = indexin(RS, ti)
  m = size(ti,1)
  alpha = ones(m)
  A = spzeros(m,m)
  logprob = zeros(maxit)
  for i in 2:maxit
    h = alpha / sum(alpha)
    Aup!(A,le_rank, re_rank, ls_rank, rs_rank, h)
    logprob[i] = summingup!(alpha, A, h)
    if abs((logprob[i]-logprob[i-1])/logprob[i-1])<reltol
      logprob = logprob[2:i]
      break
    end
  end
  return (value=ti, alpha=alpha, A=A, logprob=logprob)
end

prob2ccdf(x) = reverse(cumsum(reverse(x)))

colmarginal(x) = cumsum(vec(sum(x,dims=2)))

function randccdf(x)
  hs = [x[i] > 0 ? rand(Gamma(x[i],1)) : 0 for i in eachindex(x)]
  hs /= sum(hs)
  return prob2ccdf(hs)
end

function randcdf(x)
  hs = [x[i] > 0 ? rand(Gamma(x[i],1)) : 0 for i in eachindex(x)]
  hs /= sum(hs)
  return cumsum(hs)
end

function randprob(x)
  hs = [x[i] > 0 ? rand(Gamma(x[i],1)) : 0 for i in eachindex(x)]
  hs /= sum(hs)
  return hs
end

function randfreq(x)
  hs = [x[i] > 0 ? rand(Gamma(x[i],1)) : 0 for i in eachindex(x)]
  return cumsum(hs)
end

function sampleccdfci(x, np, p)
  seccdf = mapreduce(permutedims, vcat, [randccdf(x) for _ in 1:np])
  return  mapreduce(permutedims, vcat, [quantile(seccdf[:,i], p) for i in axes(seccdf,2)])
end

function samplefreqci(x, np, p)
  sfreq = mapreduce(permutedims, vcat, [randfreq(x) for _ in 1:np])
  return  mapreduce(permutedims, vcat, [quantile(sfreq[:,i], p) for i in axes(sfreq,2)])
end

function empiricalcdf(p,v,x)
  out = 0
  for i in eachindex(p)
    if v[i] > x
      break
    end
    out += p[i]
  end
  return out
end

function simulated_pvalue1(Rs, breaks, pos, tvalue)
  pv = Vector{Union{Missing, Float64}}(missing, length(pos))
   for i in eachindex(pos)
    fi = breaks .< pos[i]
    if any(fi)
    ind = findlast(fi)
    rs = Rs[:,ind]
    lower = ecdf(rs)(tvalue[i])
    pv[i] = 2*min(lower, 1 - lower)
    end
  end
  return pv
end

function simulated_pvalue2(Rs, breaks, pos, tvalue)
  pv = Vector{Union{Missing, Float64}}(missing, length(pos))
  #ps = collect(range(0, stop=1, length = size(Rs, 1)))
  for i in eachindex(pos)
    #fi = breaks .< pos[i]
    fi = breaks .> pos[i]
    if any(fi)
    #ind = findlast(fi)
    ind = findfirst(fi)
    rs = Rs[:,ind]
    lower = ecdf(rs)(tvalue[i])
    pv[i] = 2*min(lower, 1 - lower)
    end
  end
  return pv
end

function simulated_pvalue(Rs, tvalue)
  lower = ecdf(Rs)(tvalue)
  pv = 2*min(lower, 1 - lower)
  return pv
end

function randfreqn(α, n)
  return mapreduce(permutedims, vcat, [randfreq(α) for _ in 1:n])  
end

function takeCI(Rs, probs)
  mapreduce(permutedims, vcat,[quantile(Rs[:,i], probs) for i in axes(Rs, 2)])  
end

end