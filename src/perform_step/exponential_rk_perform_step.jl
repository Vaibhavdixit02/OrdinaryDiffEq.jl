using LinearAlgebra: axpy!

# Helper function to compute the G_nj factors for the classical ExpRK methods
@inline _compute_nl(f::SplitFunction, u, p, t, A) = f.f2(u, p, t)
@inline _compute_nl(f::ODEFunction, u, p, t, A) = f(u, p, t) - A * u
@inline _compute_nl!(G, f::SplitFunction, u, p, t, A, Au_cache) = f.f2(G, u, p, t)
@inline function _compute_nl!(G, f::ODEFunction, u, p, t, A, Au_cache)
  f(G, u, p, t)
  mul!(Au_cache, A, u)
  G .-= Au_cache
end

##########################################
# Common initializers for ExpRK integrators
function initialize!(integrator, cache::ExpRKConstantCache)
  # Pre-start fsal
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
  integrator.fsallast = zero(integrator.fsalfirst)

  # Initialize interpolation derivatives
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end
function initialize!(integrator, cache::ExpRKCache)
  # Pre-start fsal
  integrator.fsalfirst = zero(cache.rtmp)
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.fsallast = zero(integrator.fsalfirst)

  # Initialize interpolation derivatives
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

###########################################
# Classical ExpRK integrators
function perform_step!(integrator, cache::LawsonEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  nl = _compute_nl(f, uprev, p, t, A)
  @muladd v = uprev + dt * nl
  if alg.krylov
    u = expv(dt, A, v; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
  else
    exphA = cache.ops
    u = exphA * v
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::LawsonEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,G,Jcache,exphA,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : (f.jac(Jcache, uprev, p, t); Jcache) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  _compute_nl!(G, f, uprev, p, t, A, rtmp)
  @muladd @. tmp = uprev + dt*G
  if alg.krylov
    Ks, expv_cache = KsCache
    arnoldi!(Ks, f.f1, tmp; m=min(alg.m, size(f.f1,1)), norm=integrator.opts.internalnorm, 
      cache=u, iop=alg.iop)
    expv!(u,dt,Ks; cache=expv_cache)
  else
    mul!(u,exphA,tmp)
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::NorsettEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  if alg.krylov
    w = phiv(dt, A, integrator.fsalfirst, 1; m=min(alg.m, size(A,1)), 
      norm=integrator.opts.internalnorm, iop=alg.iop)
    u = uprev + dt * w[:,2]
  else
    phihA = cache.ops
    u = uprev + dt * (phihA * integrator.fsalfirst)
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::NorsettEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack rtmp,Jcache,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : (f.jac(Jcache, uprev, p, t); Jcache) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  if alg.krylov
    Ks, phiv_cache, ws = KsCache; w = ws[1]
    arnoldi!(Ks, A, integrator.fsalfirst; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, 
      cache=u, iop=alg.iop)
    phiv!(w, dt, Ks, 1; cache=phiv_cache)
    @muladd @. u = uprev + dt * @view(w[:, 2])
  else
    mul!(rtmp, cache.phihA, integrator.fsalfirst)
    @muladd @. u = uprev + dt*rtmp
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ETDRK2ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  if alg.krylov
    F1 = integrator.fsalfirst
    w1 = phiv(dt, A, F1, 2; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    U2 = uprev + dt * w1[:, 2]
    F2 = _compute_nl(f, U2, p, t + dt, A) + A * uprev
    w2 = phiv(dt, A, F2, 2; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, iop=alg.iop)
    u = uprev + dt * (w1[:, 2] - w1[:, 3] + w2[:, 3])
  else
    phi1, phi2 = cache.ops
    # The caching version uses a special formula to save computation
    G1 = f.f2(uprev, p, t)
    F1 = integrator.fsalfirst
    U2 = uprev + dt * (phi1 * F1)
    G2 = f.f2(U2, p, t + dt)
    u = U2 + dt * (phi2 * (G2 - G1))
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::ETDRK2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,F2,Jcache,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : (f.jac(Jcache, uprev, p, t); Jcache) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  if alg.krylov
    F1 = integrator.fsalfirst
    Ks, phiv_cache, ws = KsCache
    w1, w2 = ws
    # Krylov for F1
    arnoldi!(Ks, A, F1; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w1, dt, Ks, 2; cache=phiv_cache)
    # Krylov for F2
    @muladd @. tmp = uprev + dt * @view(w1[:, 2])
    _compute_nl!(F2, f, tmp, p, t + dt, A, rtmp)
    F2 .+= mul!(rtmp, A, uprev)
    arnoldi!(Ks, A, F2; m=min(alg.m, size(A,1)), norm=integrator.opts.internalnorm, cache=tmp, iop=alg.iop)
    phiv!(w2, dt, Ks, 2; cache=phiv_cache)
    # Update u
    u .= uprev
    axpy!( dt, @view(w1[:, 2]), u)
    axpy!(-dt, @view(w1[:, 3]), u)
    axpy!( dt, @view(w2[:, 3]), u)
  else
    phi1, phi2 = cache.ops
    F1 = integrator.fsalfirst
    # The caching version uses a special formula to save computation
    # Compute U2
    mul!(rtmp, phi1, F1)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U2
    # Compute G2 - G1, storing result in the cache F2
    f.f2(rtmp, uprev, p, t)
    f.f2(F2, tmp, p, t + dt)
    F2 .-= rtmp # "F2" is G2 - G1
    # Update u
    u .= tmp
    axpy!(dt, mul!(rtmp, phi2, F2), u)
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ETDRK3ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  Au = A * uprev
  F1 = integrator.fsalfirst
  if alg.krylov
    # TODO: change to named tuple in v0.7
    kwargs = [(:m, min(alg.m, size(A,1))), (:norm, integrator.opts.internalnorm), (:iop, alg.iop)]
    # Krylov on F1 (first column)
    Ks = arnoldi(A, F1; kwargs...)
    w1_half = phiv(dt/2, Ks, 1)
    w1 = phiv(dt, Ks, 3)
    U2 = uprev + dt/2 * w1_half[:, 2]
    F2 = _compute_nl(f, U2, p, t + dt/2, A) + Au
    # Krylov on F2 (second column)
    w2 = phiv(dt, A, F2, 3; kwargs...)
    U3 = uprev + dt * (2w2[:, 2] - w1[:, 2])
    F3 = _compute_nl(f, U3, p, t + dt, A) + Au
    # Krylov on F3 (third column)
    w3 = phiv(dt, A, F3, 3; kwargs...)
    u = uprev + dt * (4w1[:,4] - 3w1[:,3] + w1[:,2]
                      -8w2[:,4] + 4w2[:,3]
                      +4w3[:,4] - w3[:,3])
  else
    A21, A3, B1, B2, B3 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    U2 = uprev + dt * (A21 * F1)
    F2 = f.f2(U2, p, t + dt/2) + Au
    # stage 3
    U3 = uprev + dt * (A3 * (2F2 - F1))
    F3 = f.f2(U3, p, t + dt) + Au
    # update u
    u = uprev + dt * (B1 * F1 + B2 * F2 + B3 * F3)
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::ETDRK3Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,Au,F2,F3,Jcache,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : (f.jac(Jcache, uprev, p, t); Jcache) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  F1 = integrator.fsalfirst
  mul!(Au, A, uprev)
  halfdt = dt/2
  if alg.krylov
    Ks, phiv_cache, ws = KsCache
    w1_half, w1, w2, w3 = ws
    # TODO: change to named tuple in v0.7
    kwargs = [(:m, min(alg.m, size(A,1))), (:norm, integrator.opts.internalnorm), (:iop, alg.iop), (:cache, tmp)]
    # Krylov for F1 (first column)
    arnoldi!(Ks, A, F1; kwargs...)
    phiv!(w1_half, halfdt, Ks, 1; cache=phiv_cache)
    phiv!(w1, dt, Ks, 3; cache=phiv_cache)
    @muladd @. @views tmp = uprev + halfdt * w1_half[:, 2] # tmp is U2
    _compute_nl!(F2, f, tmp, p, t + halfdt, A, rtmp); F2 .+= Au
    # Krylov for F2 (second column)
    arnoldi!(Ks, A, F2; kwargs...)
    phiv!(w2, dt, Ks, 3; cache=phiv_cache)
    @muladd @. @views tmp = uprev + dt * (2*w2[:, 2] - w1[:, 2]) # tmp is U3
    _compute_nl!(F3, f, tmp, p, t + dt, A, rtmp); F3 .+= Au
    # Krylov for F3 (third column)
    arnoldi!(Ks, A, F3; kwargs...)
    phiv!(w3, dt, Ks, 3; cache=phiv_cache)
    # Update u
    @views @. rtmp = 4w1[:,4] - 3w1[:,3] + w1[:,2] - 8w2[:,4] + 4w2[:,3] + 4w3[:,4] - w3[:,3]
    @muladd @. u = uprev + dt * rtmp
  else
    A21, A3, B1, B2, B3 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    mul!(rtmp, A21, F1)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U2
    f.f2(F2, tmp, p, t + halfdt); F2 .+= Au
    # stage 3
    @muladd @. F3 = 2 * F2 - F1 # use F3 temporarily as cache
    mul!(rtmp, A3, F3)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U3
    f.f2(F3, tmp, p, t + dt); F3 .+= Au
    # update u
    u .= uprev
    axpy!(dt, mul!(rtmp, B1, F1), u)
    axpy!(dt, mul!(rtmp, B2, F2), u)
    axpy!(dt, mul!(rtmp, B3, F3), u)
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ETDRK4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  Au = A * uprev
  F1 = integrator.fsalfirst
  halfdt = dt/2
  if alg.krylov
    # TODO: change to named tuple in v0.7
    kwargs = [(:m, min(alg.m, size(A,1))), (:norm, integrator.opts.internalnorm), (:iop, alg.iop)]
    # Krylov on F1 (first column)
    Ks = arnoldi(A, F1; kwargs...)
    w1_half = phiv(halfdt, Ks, 1)
    w1 = phiv(dt, Ks, 3)
    U2 = uprev + halfdt * w1_half[:, 2]
    F2 = _compute_nl(f, U2, p, t + halfdt, A) + Au
    # Krylov on F2 (second column)
    Ks = arnoldi(A, F2; kwargs...)
    w2_half = phiv(halfdt, Ks, 1)
    w2 = phiv(dt, Ks, 3)
    U3 = uprev + halfdt * w2_half[:, 2]
    F3 = _compute_nl(f, U3, p, t + halfdt, A) + Au
    # Krylov on F3 (third column)
    w3 = phiv(dt, A, F3, 3; kwargs...)
    # Extra Krylov for computing F4
    rtmp = 2F3 - F1 - Au + A*U2
    wtmp = phiv(halfdt, A, rtmp, 1; kwargs...)
    U4 = U2 + halfdt * wtmp[:, 2]
    F4 = _compute_nl(f, U4, p, t + dt, A) + Au
    # Krylov on F4 (fourth column)
    w4 = phiv(dt, A, F4, 3; kwargs...)
    # update u
    u = uprev + dt * (w1[:,2] - 3w1[:,3] + 4w1[:,4] + 2w2[:,3] - 4w2[:,4] + 
                      2w3[:,3] - 4w3[:,4] + 4w4[:,4] - w4[:,3])
  else
    A21, A41, A43, B1, B2, B4 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    U2 = uprev + dt * (A21 * F1)
    F2 = f.f2(U2, p, t + halfdt) + Au
    # stage 3
    U3 = uprev + dt * (A21 * F2) # A32 = A21
    F3 = f.f2(U3, p, t + halfdt) + Au
    # stage 4
    U4 = uprev + dt * (A41 * F1 + A43 * F3)
    F4 = f.f2(U4, p, t + dt) + Au
    # update u
    u = uprev + dt * (B1 * F1 + B2 * (F2 + F3) + B4 * F4) # B3 = B2
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::ETDRK4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,Au,F2,F3,F4,Jcache,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : (f.jac(Jcache, uprev, p, t); Jcache) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  F1 = integrator.fsalfirst
  mul!(Au, A, uprev)
  halfdt = dt/2
  if alg.krylov
    Ks, phiv_cache, ws = KsCache
    w1_half, w2_half, w1, w2, w3, w4 = ws
    # TODO: change to named tuple in v0.7
    kwargs = [(:m, min(alg.m, size(A,1))), (:norm, integrator.opts.internalnorm), (:iop, alg.iop), (:cache, tmp)]
    # Krylov for F1 (first column)
    arnoldi!(Ks, A, F1; kwargs...)
    phiv!(w1_half, halfdt, Ks, 1; cache=phiv_cache)
    phiv!(w1, dt, Ks, 3; cache=phiv_cache)
    U2 = u # temporarily use u to store U2 (used in the extra Krylov step)
    @muladd @. @views U2 = uprev + halfdt * w1_half[:, 2]
    _compute_nl!(F2, f, U2, p, t + halfdt, A, rtmp); F2 .+= Au
    # Krylov for F2 (second column)
    arnoldi!(Ks, A, F2; kwargs...)
    phiv!(w2_half, halfdt, Ks, 1; cache=phiv_cache)
    phiv!(w2, dt, Ks, 3; cache=phiv_cache)
    @muladd @. @views tmp = uprev + halfdt * w2_half[:, 2] # tmp is U3
    _compute_nl!(F3, f, tmp, p, t + halfdt, A, rtmp); F3 .+= Au
    # Krylov for F3 (third column)
    arnoldi!(Ks, A, F3; kwargs...)
    phiv!(w3, dt, Ks, 3; cache=phiv_cache)
    # Extra Krylov for computing F4
    # Compute rtmp = 2F3 - F1 - Au + A*U2
    mul!(rtmp, A, U2); @. rtmp += 2F3 - F1 - Au
    arnoldi!(Ks, A, rtmp; kwargs...)
    phiv!(w1_half, halfdt, Ks, 1; cache=phiv_cache) # original w1_half is no longer needed
    @muladd @. @views tmp = U2 + halfdt * w1_half[:, 2] # tmp is U4
    _compute_nl!(F4, f, tmp, p, t + dt, A, rtmp); F4 .+= Au
    # Krylov for F4 (fourth column)
    arnoldi!(Ks, A, F4; kwargs...)
    phiv!(w4, dt, Ks, 3; cache=phiv_cache)
    # update u
    @views @. rtmp = w1[:,2] - 3w1[:,3] + 4w1[:,4] + 2w2[:,3] - 4w2[:,4] + 
                     2w3[:,3] - 4w3[:,4] + 4w4[:,4] - w4[:,3]
    @muladd @. u = uprev + dt * rtmp
  else
    A21, A41, A43, B1, B2, B4 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    mul!(rtmp, A21, F1)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U2
    f.f2(F2, tmp, p, t + halfdt); F2 .+= Au
    # stage 3
    mul!(rtmp, A21, F2) # A32 = A21
    @muladd @. tmp = uprev + dt * rtmp # tmp is U3
    f.f2(F3, tmp, p, t + halfdt); F3 .+= Au
    # stage 4
    @. tmp = uprev
    axpy!(dt, mul!(rtmp, A41, F1), tmp)
    axpy!(dt, mul!(rtmp, A43, F3), tmp) # tmp is U4
    f.f2(F4, tmp, p, t + dt); F4 .+= Au
    # update u
    u .= uprev
    axpy!(dt, mul!(rtmp, B1, F1), u)
    F2 .+= F3; axpy!(dt, mul!(rtmp, B2, F2), u) # B3 = B2
    axpy!(dt, mul!(rtmp, B4, F4), u)
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::HochOst4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = isa(f, SplitFunction) ? f.f1 : f.jac(uprev, p, t) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  Au = A * uprev
  F1 = integrator.fsalfirst
  halfdt = dt/2
  if alg.krylov
    # TODO: change to named tuple in v0.7
    kwargs = [(:m, min(alg.m, size(A,1))), (:norm, integrator.opts.internalnorm), (:iop, alg.iop)]
    # Krylov on F1 (first column)
    Ks = arnoldi(A, F1; kwargs...)
    w1_half = phiv(halfdt, Ks, 3)
    w1 =      phiv(dt,     Ks, 3)
    U2 = uprev + halfdt * w1_half[:, 2]
    F2 = _compute_nl(f, U2, p, t + halfdt, A) + Au
    # Krylov on F2 (second column)
    Ks = arnoldi(A, F2; kwargs...)
    w2_half = phiv(halfdt, Ks, 3)
    w2 =      phiv(dt,     Ks, 3)
    U3 = uprev + dt * (0.5w1_half[:,2] - w1_half[:,3] + w2_half[:,3])
    F3 = _compute_nl(f, U3, p, t + halfdt, A) + Au
    # Krylov on F3 (third column)
    Ks = arnoldi(A, F3; kwargs...)
    w3_half = phiv(halfdt, Ks, 3)
    w3 =      phiv(dt,     Ks, 3)
    U4 = uprev + dt * (w1[:,2] - 2w1[:,3] + w2[:,3] + w3[:,3])
    F4 = _compute_nl(f, U4, p, t + dt, A) + Au
    # Krylov on F4 (fourth column)
    Ks = arnoldi(A, F4; kwargs...)
    w4_half = phiv(halfdt, Ks, 3)
    w4 =      phiv(dt,     Ks, 3)
    U5 = uprev + dt * (0.5w1_half[:,2] - 0.75w1_half[:,3] + 0.5w1_half[:,4] + w1[:,4] - 0.25w1[:,3] + 
                       0.5w2_half[:,3] - w2[:,4] + 0.25w2[:,3] - 0.5w2_half[:,4] + 
                       0.5w3_half[:,3] - w2[:,4] + 0.25w3[:,3] - 0.5w3_half[:,4] + 
                       w4[:,4] - 0.25w4[:,3] - 0.25w4_half[:,3] + 0.5w4_half[:,4])
    F5 = _compute_nl(f, U5, p, t + halfdt, A) + Au
    # Krylov on F5 (fifth column)
    w5 = phiv(dt, A, F5, 3; kwargs...)
    # update u
    u = uprev + dt * (w1[:,2] - 3w1[:,3] + 4w1[:,4] - w4[:,3] + 4w4[:,4] + 4w5[:,3] - 8w5[:,4])
  else
    A21, A31, A32, A41, A42, A51, A52, A54, B1, B4, B5 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    U2 = uprev + dt * (A21 * F1)
    F2 = f.f2(U2, p, t + halfdt) + Au
    # stage 3
    U3 = uprev + dt * (A31 * F1 + A32 * F2)
    F3 = f.f2(U3, p, t + halfdt) + Au
    # stage 4
    U4 = uprev + dt * (A41 * F1 + A42 * (F2 + F3))
    F4 = f.f2(U4, p, t + dt) + Au
    # stage 5
    U5 = uprev + dt * (A51 * F1 + A52 * (F2 + F3) + A54 * F4)
    F5 = f.f2(U5, p, t + halfdt) + Au
    # update u
    u = uprev + dt * (B1 * F1 + B4 * F4 + B5 * F5)
  end

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::HochOst4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,rtmp2,Au,F2,F3,F4,F5,Jcache,KsCache = cache
  A = isa(f, SplitFunction) ? f.f1 : (f.jac(Jcache, uprev, p, t); Jcache) # get linear operator
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg

  F1 = integrator.fsalfirst
  mul!(Au, A, uprev)
  halfdt = dt/2
  if alg.krylov
    Ks, phiv_cache, ws = KsCache
    w1_half, w2_half, w3_half, w4_half, w1, w2, w3, w4, w5 = ws
    # TODO: change to named tuple in v0.7
    kwargs = [(:m, min(alg.m, size(A,1))), (:norm, integrator.opts.internalnorm), (:iop, alg.iop), (:cache, tmp)]
    # Krylov on F1 (first column)
    arnoldi!(Ks, A, F1; kwargs...)
    phiv!(w1_half, halfdt, Ks, 3; cache=phiv_cache)
    phiv!(w1,          dt, Ks, 3; cache=phiv_cache)
    @muladd @. @views tmp = uprev + halfdt * w1_half[:, 2] # tmp is U2
    _compute_nl!(F2, f, tmp, p, t + halfdt, A, rtmp); F2 .+= Au
    # Krylov on F2 (second column)
    arnoldi!(Ks, A, F2; kwargs...)
    phiv!(w2_half, halfdt, Ks, 3; cache=phiv_cache)
    phiv!(w2,          dt, Ks, 3; cache=phiv_cache)
    @muladd @. @views tmp = uprev + dt * (0.5w1_half[:,2] - w1_half[:,3] + w2_half[:,3]) # tmp is U3
    _compute_nl!(F3, f, tmp, p, t + halfdt, A, rtmp); F3 .+= Au
    # Krylov on F3 (third column)
    arnoldi!(Ks, A, F3; kwargs...)
    phiv!(w3_half, halfdt, Ks, 3; cache=phiv_cache)
    phiv!(w3,          dt, Ks, 3; cache=phiv_cache)
    @muladd @. @views tmp = uprev + dt * (w1[:,2] - 2w1[:,3] + w2[:,3] + w3[:,3]) # tmp is U4
    _compute_nl!(F4, f, tmp, p, t + dt, A, rtmp); F4 .+= Au
    # Krylov on F4 (fourth column)
    arnoldi!(Ks, A, F4; kwargs...)
    phiv!(w4_half, halfdt, Ks, 3; cache=phiv_cache)
    phiv!(w4,          dt, Ks, 3; cache=phiv_cache)
    @muladd @. @views tmp = uprev + dt * (
      0.5w1_half[:,2] - 0.75w1_half[:,3] + 0.5w1_half[:,4] + w1[:,4] - 0.25w1[:,3] + 
      0.5w2_half[:,3] - w2[:,4] + 0.25w2[:,3] - 0.5w2_half[:,4] + 
      0.5w3_half[:,3] - w2[:,4] + 0.25w3[:,3] - 0.5w3_half[:,4] + 
      w4[:,4] - 0.25w4[:,3] - 0.25w4_half[:,3] + 0.5w4_half[:,4]) # tmp is U5
    _compute_nl!(F5, f, tmp, p, t + dt, A, rtmp); F5 .+= Au
    # Krylov on F5 (fifth column)
    arnoldi!(Ks, A, F5; kwargs...)
    phiv!(w5, dt, Ks, 3; cache=phiv_cache)
    # update u
    @muladd @. @views rtmp = w1[:,2] - 3w1[:,3] + 4w1[:,4] - w4[:,3] + 4w4[:,4] + 4w5[:,3] - 8w5[:,4]
    @muladd @. u = uprev + dt * rtmp
  else
    A21, A31, A32, A41, A42, A51, A52, A54, B1, B4, B5 = cache.ops
    # stage 1 (fsaled)
    # stage 2
    mul!(rtmp, A21, F1)
    @muladd @. tmp = uprev + dt * rtmp # tmp is U2
    f.f2(F2, tmp, p, t + halfdt); F2 .+= Au
    # stage 3
    mul!(rtmp, A31, F1); mul!(rtmp2, A32, F2); rtmp .+= rtmp2
    @muladd @. tmp = uprev + dt * rtmp # tmp is U3
    f.f2(F3, tmp, p, t + halfdt); F3 .+= Au
    # stage 4
    F2 .+= F3 # F2 now stores F2 + F3
    mul!(rtmp, A41, F1); mul!(rtmp2, A42, F2); rtmp .+= rtmp2
    @muladd @. tmp = uprev + dt * rtmp # tmp is U4
    f.f2(F4, tmp, p, t + dt); F4 .+= Au
    # stage 5
    mul!(rtmp, A51, F1)
    mul!(rtmp2, A52, F2); rtmp .+= rtmp2
    mul!(rtmp2, A54, F4); rtmp .+= rtmp2
    @muladd @. tmp = uprev + dt * rtmp # tmp is U5
    f.f2(F5, tmp, p, t + halfdt); F5 .+= Au
    # update u
    mul!(rtmp, B1, F1)
    mul!(rtmp2, B4, F4); rtmp .+= rtmp2
    mul!(rtmp2, B5, F5); rtmp .+= rtmp2
    @muladd @. u = uprev + dt * rtmp
  end

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

#############################################
# EPIRK integrators
function perform_step!(integrator, cache::Exp4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = f.jac(uprev, p, t)
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  f0 = integrator.fsalfirst # f(uprev) is fsaled
  ts = [dt/3, 2dt/3, dt]
  kwargs = [(:tol, integrator.opts.reltol), (:iop, alg.iop), (:norm, integrator.opts.internalnorm), (:adaptive, true)]

  # Krylov for f(uprev)
  B1 = [zero(f0) f0]
  K1 = phiv_timestep(ts, A, B1; kwargs...) # tϕ(tA)f0
  @inbounds for i = 1:3
    K1[:,i] ./= ts[i]
  end
  w4 = K1 * [-7/300, 97/150, -37/300]
  u4 = uprev + dt * w4
  d4 = f(u4, p, t+dt) - f0 - dt * (A*w4) # TODO: what should be the time?
  # Krylov for the first remainder d4
  B2 = [zero(d4) d4]
  K2 = phiv_timestep(ts, A, B2; kwargs...)
  @inbounds for i = 1:3
    K2[:,i] ./= ts[i]
  end
  w7 = K1 * [59/300, -7/75, 269/300] + K2 * [2/3, 2/3, 2/3]
  u7 = uprev + dt * w7
  d7 = f(u7, p, t+dt) - f0 - dt * (A*w7)
  # Krylov for the second remainder d7
  B3 = [zero(d7) d7]
  k7 = phiv_timestep(ts[1], A, B3; kwargs...)
  k7 ./= ts[1]
  # Update u
  u = uprev + dt * (K1[:,3] + K2[:,1] - 4/3*K2[:,2] + K2[:,3] + k7/6)

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::Exp4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,rtmp2,K,A,B,KsCache = cache
  f.jac(A, uprev, p, t)
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  f0 = integrator.fsalfirst # f(u0) is fsaled
  ts = [dt/3, 2dt/3, dt]
  kwargs = [(:tol, integrator.opts.reltol), (:iop, alg.iop), (:norm, integrator.opts.internalnorm), 
    (:adaptive, true), (:caches, KsCache)]

  # Krylov for f(uprev)
  B[:, 2] .= f0
  phiv_timestep!(K, ts, A, B; kwargs...)
  @inbounds for i = 1:3
    K[:,i] ./= ts[i]
  end
  mul!(rtmp, K, [-7/300, 97/150, -37/300]) # rtmp is now w4
  @muladd @. tmp = uprev + dt * rtmp # tmp is now u4
  mul!(rtmp2, A, rtmp)
  f(rtmp, tmp, p, t+dt) # TODO: what should be the time?
  @muladd @. @view(B[:,2]) = rtmp - f0 - dt * rtmp2 # B[:,2] is now d4
  # Partially update entities that use k1, k2, k3
  mul!(rtmp, K, [59/300, -7/75, 269/300]) # rtmp is now w7
  @muladd @. u = uprev + dt * @view(K[:,3])
  # Krylov for the first remainder d4
  phiv_timestep!(K, ts, A, B; kwargs...)
  @inbounds for i = 1:3
    K[:,i] ./= ts[i]
  end
  mul!(rtmp2, K, [2/3, 2/3, 2/3]); rtmp .+= rtmp2 # w7 fully updated
  @muladd @. tmp = uprev + dt * rtmp # tmp is now u7
  mul!(rtmp2, A, rtmp)
  f(rtmp, tmp, p, t+dt) # TODO: what should be the time?
  @muladd @. @view(B[:,2]) = rtmp - f0 - dt * rtmp2 # B[:,2] is now d7
  # Partially update entities that use k4, k5, k6
  mul!(rtmp, K, [1.0, -4/3, 1.0])
  axpy!(dt, rtmp, u)
  # Krylov for the second remainder d7
  k7 = @view(K[:, 1])
  phiv_timestep!(k7, ts[1], A, B; kwargs...)
  k7 ./= ts[1]
  axpy!(dt/6, k7, u)

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::EPIRK4s3AConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = f.jac(uprev, p, t)
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  f0 = integrator.fsalfirst # f(uprev) is fsaled
  kwargs = [(:tol, integrator.opts.reltol), (:iop, alg.iop), (:norm, integrator.opts.internalnorm), (:adaptive, true)]

  # Compute U2 and U3 vertically
  K = phiv_timestep([dt/2, 2dt/3], A, [zero(f0) f0]; kwargs...)
  U2 = uprev + K[:, 1] 
  U3 = uprev + K[:, 2]
  R2 = f(U2, p, t + dt/2)  - f0 - A*K[:, 1] # remainder of U2
  R3 = f(U3, p, t + 2dt/3) - f0 - A*K[:, 2] # remainder of U3
  
  # Update u (horizontally)
  B = zero(eltype(f0), length(f0), 5)
  B[:, 2] = f0
  B[:, 4] = (32R2 - 13.5R3) / dt^2
  B[:, 5] = (-144R2 + 81R3) / dt^3
  u = uprev + phiv_timestep(dt, A, B; kwargs...)

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::EPIRK4s3ACache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,rtmp2,K,A,B,KsCache = cache
  f.jac(A, uprev, p, t)
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  f0 = integrator.fsalfirst # f(u0) is fsaled
  kwargs = [(:tol, integrator.opts.reltol), (:iop, alg.iop), (:norm, integrator.opts.internalnorm), 
    (:adaptive, true), (:caches, KsCache)]

  # Compute U2 and U3 vertically
  B[:, 2] .= f0
  phiv_timestep!(K, [dt/2, 2dt/3], A, @view(B[:, 1:2]); kwargs...)
  ## U2 and R2
  @. tmp = uprev + @view(K[:, 1]) # tmp is now U2
  f(rtmp, tmp, p, t + dt/2); mul!(rtmp2, A, @view(K[:, 1]))
  @. rtmp = rtmp - f0 - rtmp2 # rtmp is now R2
  B[:, 4] .= (32/dt^2) * rtmp
  B[:, 5] .= (-144/dt^3) * rtmp
  ## U3 and R3
  @. tmp = uprev + @view(K[:, 2]) # tmp is now U3
  f(rtmp, tmp, p, t + 2dt/3); mul!(rtmp2, A, @view(K[:, 2]))
  @. rtmp = rtmp - f0 - rtmp2 # rtmp is now R3
  B[:, 4] .-= (13.5/dt^2) * rtmp
  B[:, 5] .+= (81/dt^3) * rtmp
  
  # Update u
  du = @view(K[:, 1])
  phiv_timestep!(du, dt, A, B; kwargs...)
  @. u = uprev + du

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::EPIRK4s3BConstantCache, repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  A = f.jac(uprev, p, t)
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  f0 = integrator.fsalfirst # f(uprev) is fsaled
  kwargs = [(:tol, integrator.opts.reltol), (:iop, alg.iop), (:norm, integrator.opts.internalnorm), (:adaptive, true)]

  # Compute U2 and U3 vertically
  K = phiv_timestep([dt/2, 3dt/4], A, [zero(f0) zero(f0) f0]; kwargs...)
  K[:, 1] .*= 8 / (3*dt)
  K[:, 2] .*= 16 / (9*dt)
  U2 = uprev + K[:, 1]
  U3 = uprev + K[:, 2]
  R2 = f(U2, p, t + dt/2)  - f0 - A*K[:, 1] # remainder of U2
  R3 = f(U3, p, t + 3dt/4) - f0 - A*K[:, 2] # remainder of U3
  
  # Update u (horizontally)
  B = zero(eltype(f0), length(f0), 5)
  B[:, 2] = f0
  B[:, 4] = (54R2 - 16R3) / dt^2
  B[:, 5] = (-324R2 + 144R3) / dt^3
  u = uprev + phiv_timestep(dt, A, B; kwargs...)

  # Update integrator state
  integrator.fsallast = f(u, p, t + dt)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function perform_step!(integrator, cache::EPIRK4s3BCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,rtmp,rtmp2,K,A,B,KsCache = cache
  f.jac(A, uprev, p, t)
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  f0 = integrator.fsalfirst # f(u0) is fsaled
  kwargs = [(:tol, integrator.opts.reltol), (:iop, alg.iop), (:norm, integrator.opts.internalnorm), 
    (:adaptive, true), (:caches, KsCache)]

  # Compute U2 and U3 vertically
  fill!(@view(B[:, 2]), zero(eltype(B)))
  B[:, 3] .= f0
  phiv_timestep!(K, [dt/2, 3dt/4], A, @view(B[:, 1:3]); kwargs...)
  K[:, 1] .*= 8 / (3*dt)
  K[:, 2] .*= 16 / (9*dt)
  ## U2 and R2
  @. tmp = uprev + @view(K[:, 1]) # tmp is now U2
  f(rtmp, tmp, p, t + dt/2); mul!(rtmp2, A, @view(K[:, 1]))
  @. rtmp = rtmp - f0 - rtmp2 # rtmp is now R2
  B[:, 4] .= (54/dt^2) * rtmp
  B[:, 5] .= (-324/dt^3) * rtmp
  ## U3 and R3
  @. tmp = uprev + @view(K[:, 2]) # tmp is now U3
  f(rtmp, tmp, p, t + 3dt/4); mul!(rtmp2, A, @view(K[:, 2]))
  @. rtmp = rtmp - f0 - rtmp2 # rtmp is now R3
  B[:, 4] .-= (16/dt^2) * rtmp
  B[:, 5] .+= (144/dt^3) * rtmp
  
  # Update u
  fill!(@view(B[:, 3]), zero(eltype(B)))
  B[:, 2] .= f0
  du = @view(K[:, 1])
  phiv_timestep!(du, dt, A, B; kwargs...)
  @. u = uprev + du

  # Update integrator state
  f(integrator.fsallast, u, p, t + dt)
  # integrator.k is automatically set due to aliasing
end

######################################################
# Multistep exponential integrators
function initialize!(integrator,cache::ETD2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Pre-start fsal
  lin = integrator.f.f1(integrator.uprev,integrator.p,integrator.t)
  nl = integrator.f.f2(integrator.uprev,integrator.p,integrator.t)
  nlprev = zero(nl) # to be computed in the first iteration via ETD1
  integrator.fsalfirst = ETD2Fsal(lin, nl, nlprev)
    
  # Avoid undefined entries if k is an array of arrays
  rate_prototype = lin
  integrator.fsallast = ETD2Fsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator,cache::ETD2ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,f,p = integrator
  @unpack lin,nl,nlprev = integrator.fsalfirst
  @unpack exphA,phihA,B1,B0 = cache
  integrator.k[1] = lin + nl

  if integrator.iter == 1 # ETD1 for initial step
    u = exphA*uprev + dt*(phihA*nl)
  else
    u = exphA*uprev + dt*(B1*nl + B0*nlprev)
  end
  integrator.u = u

  # Push the fsal at t+dt
  nlprev = nl
  lin = f.f1(u,p,t+dt)
  nl = f.f2(u,p,t+dt)
  integrator.k[2] = lin + nl
  @pack integrator.fsallast = lin, nl, nlprev
end

function initialize!(integrator, cache::ETD2Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  rate_prototype = cache.rtmp1
  
  # Pre-start fsal
  integrator.fsalfirst = ETD2Fsal(rate_prototype)
  @unpack lin,nl = integrator.fsalfirst
  integrator.f.f1(lin,integrator.uprev,integrator.p,integrator.t)
  integrator.f.f2(nl,integrator.uprev,integrator.p,integrator.t)
    
  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = ETD2Fsal(rate_prototype)
  integrator.k[1] = lin + nl
  integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::ETD2Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack lin,nl,nlprev = integrator.fsalfirst
  @unpack utmp,rtmp1,rtmp2,exphA,phihA,B1,B0 = cache
  @. integrator.k[1] = lin + nl

  if integrator.iter == 1 # ETD1 for initial step
    mul!(utmp, exphA, uprev)
    mul!(rtmp1, phihA, nl)
    @muladd @. u = utmp + dt*rtmp1
  else
    mul!(utmp, exphA, uprev)
    mul!(rtmp1, B1, nl)
    mul!(rtmp2, B0, nlprev)
    @muladd @. u = utmp + dt*(rtmp1 + rtmp2)
  end

  # Push the fsal at t+dt
  fsallast = integrator.fsallast
  fsallast.nlprev .= nl
  f.f1(fsallast.lin,u,p,t+dt)
  f.f2(fsallast.nl,u,p,t+dt)
  @. integrator.k[2] = fsallast.lin + fsallast.nl
end
