function potential_hartree_g!(g::Vector{Float64},
                              rhog::Vector{Complex{Float64}},
                              vhg::Vector{Complex{Float64}})
    vhg[1] = 0.0 + 0.0im
    ng = size(g)[1]
    for i in 2:ng
        vhg[i] = 4π * rhog[i] / g[i] / g[i]
    end
end
