const soh     = Libdl.dlopen(Pkg.dir("IScatterSpectrum")*"/libiscatspe")
const plascb  = Libdl.dlsym(soh, :plas_)
const plasmaf = Libdl.dlsym(soh, :plasma_)
const pf32    = cglobal(plascb, Float32)
const pi32    = cglobal(plascb, Int32)

zr     = unsafe_wrap(Vector{Float32}, pf32, 100)
zi     = unsafe_wrap(Vector{Float32}, pf32+4*100, 100)
x      = unsafe_wrap(Vector{Float32}, pf32+4*200, 100)
y      = unsafe_wrap(Vector{Float32}, pf32+4*300, 1)
f      = unsafe_wrap(Vector{Float32}, pf32+4*301, 120)
scalef = unsafe_wrap(Vector{Float32}, pf32+4*421, 1)
nx     = unsafe_wrap(Vector{Int32}, pi32+4*422, 1)

function plasma!(freq::Vector{Float32})
    for (k, f) in enumerate(freq)
        unsafe_store!(pf32, f, 301+k)
    end
    unsafe_store!(pi32, Int32(length(freq)), 423)
    return
end

function plasma!(scalef::Float32=1f0, y::Float32=0f0)
    unsafe_store!(pf32,      y, 301)
    unsafe_store!(pf32, scalef, 422)
    ccall(plasmaf, Void, ())
    nx = unsafe_load(pi32, 423)
    hcat([unsafe_load(pf32, k) for k in 1:nx],
         [unsafe_load(pf32, k+100) for k in 1:nx])
end
