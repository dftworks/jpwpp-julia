include("lattice.jl")
include("fftgrid.jl")
include("gvector.jl")
include("pwdensity.jl")
include("kgrid.jl")
include("pwbasis.jl")
include("atom.jl")
include("hartree.jl")
include("atompspfhi.jl")
include("initrho.jl")

const BOHR_TO_ANG = 0.529177249
const ANG_TO_BOHR = 1 / BOHR_TO_ANG

ecut_rho = 4.0 # 40 Hartree
ecut_wfc = 1.0 # 10 Hartree

latt_const = 5.43 * ANG_TO_BOHR

latt = Lattice( latt_const * [0.5, 0.5, 0.0],
                latt_const * [0.0, 0.5, 0.5],
                latt_const * [0.5, 0.0, 0.5] )

println("\nreal space\n", latt)
println("unit cell volume = ", volume(latt))

blatt = reciprocal_lattice(latt)

println("\nreciprocal space\n", blatt)
println("Brillouin zone volume = ", volume(blatt))

println("\nverify the lattices in direct and reciprocal spaces")
println("latt_volume * blatt_volume = ", volume(latt) * volume(blatt))
println("(2*pi)^3                   = ", (2*pi)^3)

fftgrid = FFTGrid(latt, ecut_rho) # 40 Hartree

println("\n", fftgrid)

gvec = GVector(fftgrid, blatt)

println("\n", gvec)

pwden = PWDensity(ecut_rho, gvec)

println("\nNo of G for density expansion : ", pwden.npw)
println("No of G shell : ", pwden.nshell)

kgrid = KGrid(2, 1, 1)

pwbasis_array = Vector{PWBasis}(undef, kgrid.nkpt)

for ik in 1:kgrid.nkpt
    fk = kgrid.kpts[:,ik]
    ck = frac_to_cart(fk, blatt)
    pwbasis_array[ik] = PWBasis(ck, ecut_wfc, gvec)
end

for i in 1:kgrid.nkpt
    println("\nk - ", i)
    println("Fractional : ", kgrid.kpts[:,i])
    println("Cartesion  : ", frac_to_cart(kgrid.kpts[:,i], blatt))
    println("No of PWs  : ", pwbasis_array[i].npw)
end

atoms = Vector{Atom}(undef, 2)
atoms[1] = Atom("Si", 0.0, 0.0, 0.0)
atoms[2] = Atom("Si", 0.25, 0.25, 0.25)

println("\nAtoms:")
for at in atoms
    println(at)
end

atompsp = AtomPSPFHI("14-Si.LDA.fhi")

rho_gshell_1d = Vector{Complex{Float64}}(undef, pwden.nshell)

ρ_gshell_1d = compute_rho_of_gshell(atompsp.rhops, atompsp.rad, pwden.gshell, volume(latt))

vh_gshell_1d = Vector{Complex{Float64}}(undef, pwden.nshell)

potential_hartree_g!(pwden.gshell, ρ_gshell_1d, vh_gshell_1d)

sfact = compute_structure_factor(gvec.miller, atoms)
println(sfact)
