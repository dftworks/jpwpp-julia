struct KGrid
    kpts :: Matrix{Float64}
    nkpt :: Int
end
function KGrid(nk1::Int, nk2::Int, nk3::Int) ::KGrid
    nkpt = nk1 * nk2 * nk3
    kpts = zeros(3,nkpt)
    n = 0
    for ik1 in 1:nk1
        for ik2 in 1:nk2
            for ik3 in 1:nk3
                n += 1
                k1 = monkhorst_pack(ik1, nk1)
                k2 = monkhorst_pack(ik2, nk2)
                k3 = monkhorst_pack(ik3, nk3)
                kpts[:,n] = [k1, k2, k3]
            end
        end
    end
    KGrid(kpts, nkpt)
end
function monkhorst_pack(r, q)
    (2r - q - 1) / 2q
end
function frac_to_cart(fk::Vector{Float64}, blatt::Lattice) ::Vector{Float64}
    ck1 = fk[1] * blatt.v1[1] + fk[2] * blatt.v2[1] + fk[3] * blatt.v3[1]
    ck2 = fk[1] * blatt.v1[2] + fk[2] * blatt.v2[2] + fk[3] * blatt.v3[2]
    ck3 = fk[1] * blatt.v1[3] + fk[2] * blatt.v2[3] + fk[3] * blatt.v3[3]
    [ck1, ck2, ck3]
end
