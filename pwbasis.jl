include("gvector.jl")

struct PWBasis
    npw :: Int
    gindex :: Vector{Int}
    kg :: Vector{Float64}
end
function PWBasis(ck::Vector{Float64}, ecut::Float64, gvec::GVector) ::PWBasis
    npw = get_n_plane_waves(gvec, ecut, ck)
    gindex = get_plane_waves(gvec, ecut, ck, npw)
    kg = zeros(npw)
    for i in 1:npw
        kg[i] = norm(ck+gvec.cart[:,gindex[i]])
    end
    PWBasis(npw, gindex, kg)
end
