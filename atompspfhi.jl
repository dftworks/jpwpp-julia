using Printf
include("simpson_log.jl")

struct AtomPSPFHI
    zatom :: Float64
    zion :: Float64
    pspcod :: Int
    pspxc :: Int
    lmax :: Int
    lloc :: Int
    mmax :: Int
    amesh :: Float64
    rad :: Vector{Float64}
    wfc :: Matrix{Float64}
    pot :: Matrix{Float64}
    # Kleinman-Bylander nonlocal parts
    lbeta :: Vector{Int}
    beta :: Matrix{Float64}
    dfact :: Vector{Float64}
    rhops :: Vector{Float64}
end
function AtomPSPFHI(pspfile :: String) ::AtomPSPFHI
    lines = readlines(pspfile)

    zatom, zion = parse.(Float64, split(lines[2])[1:2])
    pspcod, pspxc, lmax, lloc, mmax = parse.(Int, split(lines[3])[1:5])

    wfc = zeros(Float64, mmax, lmax+1)
    pot = zeros(Float64, mmax, lmax+1)
    rad = zeros(Float64, mmax)

    amesh = parse.(Float64, split(lines[19])[2])
    
    n1 = 19
    for l in 0:lmax
        n2 = n1 + mmax
        for i in 1:mmax
            rad[i], wfc[i,l+1], pot[i,l+1] =
                parse.(Float64, split(lines[n1+i])[2:4])
        end
        n1 = n2 + 1
    end

    #
    lbeta = zeros(Int, lmax)
    n = 1
    for l in 0:lmax
        if l != lloc
            lbeta[n] = l
            n += 1
        end
    end
    
    beta = zeros(Float64, mmax, lmax)
    dfact = zeros(Float64, lmax)
    for i in 1:lmax
        l = lbeta[i]
        beta[:,i] = (pot[:,l+1] - pot[:,lloc+1]) .* wfc[:,l+1]
        work = beta[:,i] .* wfc[:,l+1]
        dfact[i] = 1.0 / simpson_log(work, rad)
    end

    rhops = zeros(mmax)
    rhops[:] += wfc[:,1] .* wfc[:,1] * 2.0 # l = 0, s orbital
    rhops[:] += wfc[:,2] .* wfc[:,2] * 2.0 # l = 1, p orbital
    
    AtomPSPFHI(zatom, zion, pspcod, pspxc, lmax, lloc, mmax,
               amesh, rad, wfc, pot, lbeta, beta, dfact, rhops)
end
