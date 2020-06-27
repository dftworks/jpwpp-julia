include("gvector.jl")

struct PWDensity
    npw :: Int32
    gindex :: Vector{Int32}
    g :: Vector{Float64}
    nshell :: Int32
    gshell :: Vector{Float64}
    gshell_index :: Vector{Int32}
end
function PWDensity(ecut::Float64, gvec::GVector) ::PWDensity
    npw = get_n_plane_waves(gvec, ecut, [0.0, 0.0, 0.0]) 
    gindex = get_plane_waves(gvec, ecut, [0.0, 0.0, 0.0], npw)
    g = zeros(npw)
    for i in 1:npw
        g[i] = norm(gvec.cart[:,gindex[i]])
    end
    
    nshell = get_n_g_shell(g)
    gshell, gshell_index = get_g_shell_and_index(g, nshell, npw)
    PWDensity(npw, gindex, g, nshell, gshell, gshell_index)
end
function get_n_g_shell(g::Vector{Float64}) ::Int32
    ng = 0
    gnorm = -1.0
    for x in g
        if abs(x - gnorm) > 1E-10
            gnorm = x
            ng += 1
        end
    end
    ng
end
function get_g_shell_and_index(g::Vector{Float64}, nshell::Int32, npw::Int32)
    #::(Vector{Float64},Vector{Int32})
    ng = 0
    gnorm = -1.0

    shell = zeros(nshell)
    index = zeros(npw)

    for i in 1:size(g)[1]
        x = g[i]
        if abs(x - gnorm) > 1E-10
            gnorm = x
            ng += 1
            shell[ng] = gnorm
        end
        index[i] = ng
    end
    shell, index
end
