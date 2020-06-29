include("simpson_log.jl")

function compute_rho_of_gshell(rhops::Vector{Float64}, rad::Vector{Float64},
                          gshell::Vector{Float64}, volume::Float64) ::Vector{Complex{Float64}}
    nshell = size(gshell)[1]

    rhog = zeros(nshell)

    mmax = size(rhops)[1]
    work = zeros(mmax)

    rhog[1] = simpson_log(rhops, rad)
    
    for iw in 2:nshell
        for i in 1:mmax
            gr = gshell[iw] * rad[i]
            work[i] = rhops[i] * sin(gr) / gr
        end
        rhog[iw] = simpson_log(work, rad)
    end
    
    rhog /= volume

    rhog
end

function compute_structure_factor(miller::Matrix{Int32}, atoms::Vector{Atom})::Vector{Complex{Float64}}
    ng = size(miller)[2]
    s :: Vector{Complex{Float64}} = zeros(ng)
    for i in 1:ng
        g = miller[:,i]
        for atom in atoms
            gr = g[1] * atom.x + g[2] * atom.y + g[3] * atom.z
            s[i] += exp(-im*2*pi*gr)
        end
    end
    s
end
