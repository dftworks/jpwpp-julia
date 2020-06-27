using Printf
using LinearAlgebra

struct Lattice
    v1 :: Vector{Float64}
    v2 :: Vector{Float64}
    v3 :: Vector{Float64}
end
function Base.show(io::IO, lattice::Lattice)
    println("Lattice:")
    @printf "%20.12f%20.12f%20.12f\n" lattice.v1[1] lattice.v1[2] lattice.v1[3]
    @printf "%20.12f%20.12f%20.12f\n" lattice.v2[1] lattice.v2[2] lattice.v2[3]
    @printf "%20.12f%20.12f%20.12f" lattice.v3[1] lattice.v3[2] lattice.v3[3]
end
function volume(latt::Lattice) ::Float64
    dot(latt.v1, cross(latt.v2, latt.v3))
end
function reciprocal_lattice(latt::Lattice) ::Lattice
    factor = 2π / volume(latt)
    b1 = cross(latt.v2, latt.v3) * factor
    b2 = cross(latt.v3, latt.v1) * factor
    b3 = cross(latt.v1, latt.v2) * factor
    Lattice(b1, b2, b3)
end
