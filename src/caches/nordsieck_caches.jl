# TODO: Optimize cache size
mutable struct AN5ConstantCache{zType,lType,dtsType,dType,tsit5Type} <: OrdinaryDiffEqConstantCache
  # `z` is the Nordsieck vector
  z::zType
  # `l` is used for the corrector iteration
  l::Vector{lType}
  # `m` is a tmp vector that is used for calculating `l`
  m::Vector{lType}
  # `c_LTE` is used for the error estimation for the current order
  c_LTE::lType
  c_conv::lType
  # `dts` stores `dt`s
  dts::dtsType
  # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
  Δ::dType
  # `Tsit5` for the first step
  tsit5tab::tsit5Type
  order::Int
end

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  N = 5
  z = [zero(rate_prototype) for i in 1:N+1]
  Δ = u
  l = zeros(tTypeNoUnits,N+1); m = zeros(l)
  c_LTE = c_conv = zero(tTypeNoUnits)
  dts = zeros(typeof(dt), 6)
  tsit5tab = Tsit5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  AN5ConstantCache(z,l,m,c_LTE,c_conv,dts,Δ,tsit5tab,1)
end

mutable struct AN5Cache{uType,dType,rateType,zType,lType,dtsType,tsit5Type} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  Δ::dType
  # Error estimation
  atmp::dType
  fsalfirst::rateType
  ratetmp::rateType
  # `z` is the Nordsieck vector
  z::zType
  # `l` is used for the corrector iteration
  l::Vector{lType}
  # `m` is a tmp vector that is used for calculating `l`
  m::Vector{lType}
  # `c_LTE` is used for the error estimation for the current order
  c_LTE::lType
  c_conv::lType
  # `dts` stores `dt`s
  dts::dtsType
  # `Tsit5` for the first step
  tsit5cache::tsit5Type
  order::Int
end

u_cache(c::AN5Cache) = ()
du_cache(c::AN5Cache) = (c.fsalfirst,)

function alg_cache(alg::AN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  #################################################
  # Tsit5
  tab = Tsit5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  # Cannot alias pointers, since we have to use `k`s to start the Nordsieck vector
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u)); tmp = similar(u)
  tsit5cache = Tsit5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
  #################################################
  N = 5
  Δ = similar(atmp)
  l = zeros(tTypeNoUnits,N+1); m = zeros(l)
  c_LTE = c_conv = zero(tTypeNoUnits)
  dts = zeros(typeof(dt), 6)
  fsalfirst = zeros(rate_prototype)
  z = [zeros(rate_prototype) for i in 1:N+1]
  for i in 1:N+1
    z[i] = zeros(rate_prototype)
  end
  ratetmp = zeros(rate_prototype)

  AN5Cache(u,uprev,tmp,Δ,atmp,fsalfirst,ratetmp,
           z,l,m,c_LTE,c_conv,dts,
           tsit5cache, 1)
end

mutable struct JVODEConstantCache{zType,lType,dtType,dType,tsit5Type,etaType} <: OrdinaryDiffEqConstantCache
  # `z` is the Nordsieck vector
  z::zType
  # `l` is used for the corrector iteration
  l::Vector{lType}
  # `m` is a tmp vector that is used for calculating `l`
  m::Vector{lType}
  # `c_LTE₊₁` is used for the error estimation for the current order + 1
  c_LTE₊₁::lType
  # `c_LTE` is used for the error estimation for the current order
  c_LTE::lType
  # `c_LTE₋₁` is used for the error estimation for the current order - 1
  c_LTE₋₁::lType
  # `c_conv` is used in convergence test
  c_conv::lType
  # `c_𝒟` is used to get the order q+2 derivative vector
  c_𝒟::lType
  prev_𝒟::lType
  # `dts` stores `dt`s
  dts::Vector{dtType}
  # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
  Δ::dType
  # `Tsit5` for the first step
  tsit5tab::tsit5Type
  # `z[Δind]` is scaled `Δ` vector
  Δind::Int
  L::Int
  # same with `order` or `q`
  order::Int
  nextorder::Int
  # number of steps to take before considering to change order
  n_wait::Int
  # `η` is `dtₙ₊₁/dtₙ`
  η  ::etaType
  ηq ::etaType
  η₊₁::etaType
  η₋₁::etaType
  maxη::etaType
  dtscale::dtType
end

function alg_cache(alg::JVODE,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  if alg.algorithm == :Adams
    N = 12
  elseif alg.algorithm == :BDF
    N = 6
  else
    error("Algorithm must be :BDF or :Adams")
  end
  z = [rate_prototype for i in 1:N+1]
  Δ = u
  l = zeros(tTypeNoUnits, N+1); m = zeros(l)
  c_LTE₊₁ = c_LTE = c_LTE₋₁ = c_conv = c_𝒟 = prev_𝒟 = zero(tTypeNoUnits)
  dts = zeros(typeof(dt),N+1)
  tsit5tab = Tsit5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  η = zero(dt/dt)
  JVODEConstantCache(z,l,m,
                     c_LTE₊₁,c_LTE,c_LTE₋₁,c_conv,c_𝒟 ,prev_𝒟,
                     dts,Δ,tsit5tab,2,1,1,2,η,η,η,η,η,dt)
end

mutable struct JVODECache{uType,rateType,zType,lType,dtType,dType,etaType,tsit5Type} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  fsalfirst::rateType
  ratetmp::rateType
  # `z` is the Nordsieck vector
  z::zType
  # `l` is used for the corrector iteration
  l::Vector{lType}
  # `m` is a tmp vector that is used for calculating `l`
  m::Vector{lType}
  # `c_LTE₊₁` is used for the error estimation for the current order + 1
  c_LTE₊₁::lType
  # `c_LTE` is used for the error estimation for the current order
  c_LTE::lType
  # `c_LTE₋₁` is used for the error estimation for the current order - 1
  c_LTE₋₁::lType
  # `c_conv` is used in convergence test
  c_conv::lType
  # `c_𝒟` is used to get the order q+2 derivative vector
  c_𝒟::lType
  prev_𝒟::lType
  # `dts` stores `dt`s
  dts::Vector{dtType}
  # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
  Δ::dType
  # Error estimation
  atmp::dType
  # `Tsit5` for the first step
  tsit5cache::tsit5Type
  L::Int
  # same with `order` or `q`
  order::Int
  nextorder::Int
  # number of steps to take before considering to change order
  n_wait::Int
  # `η` is `dtₙ₊₁/dtₙ`
  η  ::etaType
  ηq ::etaType
  η₊₁::etaType
  η₋₁::etaType
  maxη::etaType
  dtscale::dtType
end

u_cache(c::JVODECache) = ()
du_cache(c::JVODECache) = (c.fsalfirst,)

function alg_cache(alg::JVODE,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  if alg.algorithm == :Adams
    N = 12
  elseif alg.algorithm == :BDF
    N = 6
  else
    error("Algorithm must be :BDF or :Adams")
  end
  #################################################
  # Tsit5
  # Cannot alias pointers, since we have to use `k`s to start the Nordsieck vector
  tab = Tsit5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u)); tmp = similar(u)
  tsit5cache = Tsit5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
  #################################################
  fsalfirst = zeros(rate_prototype)
  Δ = similar(u,uEltypeNoUnits,indices(u))
  l = zeros(tTypeNoUnits, N+1); m = zeros(l)
  c_LTE₊₁ = c_LTE = c_LTE₋₁ = c_conv = c_𝒟 = prev_𝒟 = zero(tTypeNoUnits)
  dts = zeros(typeof(dt),N+1)
  η = zero(dt/dt)
  #################################################
  # Nordsieck Vector
  z = [zeros(rate_prototype) for i in 1:N+1]
  #################################################
  JVODECache(u,uprev,tmp,fsalfirst,ratetmp,
             z, l, m,
             c_LTE₊₁, c_LTE, c_LTE₋₁, c_conv, c_𝒟 , prev_𝒟,
             dts, Δ, atmp, tsit5cache, 2, 1, 1, 2, η, η, η, η, η, dt)
end
