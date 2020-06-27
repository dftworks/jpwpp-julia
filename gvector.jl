using LinearAlgebra

#include("fftgrid.jl")

struct GVector
    miller :: Matrix{Int32}       # {3,ng}
    cart :: Matrix{Float64}     # {3,ng}
end
function GVector(fftgrid::FFTGrid, blatt::Lattice)::GVector
    miller = make_miller(fftgrid)
    cart = make_cart(miller, blatt)

    # get the permutation that sorts the gnorm
    # the same permutation is used to sort the miller and cart
    ng = size(cart)[2]
    gnorm = zeros(ng)
    for i in 1:ng
        gnorm[i] = norm(cart[:,i])
    end
    p = sortperm(gnorm)

    miller = miller[:,p]
    cart = cart[:,p]
    
    GVector(miller, cart)
end
function Base.show(io::IO, gvec::GVector)
    print(io, "No of Gvectors: ", size(gvec.cart)[2])
end
function make_miller(fftgrid::FFTGrid) ::Matrix{Int32}
    nfft = fftgrid.n1 * fftgrid.n2 * fftgrid.n3

    i1 = fft_left_end(fftgrid.n1)
    i2 = fft_right_end(fftgrid.n1)

    j1 = fft_left_end(fftgrid.n2)
    j2 = fft_right_end(fftgrid.n2)

    k1 = fft_left_end(fftgrid.n3)
    k2 = fft_right_end(fftgrid.n3)

    miller = zeros(3,nfft)

    ig = 1

    for i in i1:i2
        for j in j1:j2
            for k in k1:k2
                miller[:,ig] = [i j k]
                ig += 1
            end
        end
    end
    miller
end
function make_cart(miller::Matrix{Int32}, blatt::Lattice) ::Matrix{Float64}
    nfft = size(miller)[2]

    b1 = blatt.v1
    b2 = blatt.v2
    b3 = blatt.v3

    cart = zeros(3,nfft)
    for i in 1:nfft 
        n1 = miller[1,i]
        n2 = miller[2,i]
        n3 = miller[3,i]

        cart[1,i] = n1 * b1[1] + n2 * b2[1] + n3 * b3[1]
        cart[2,i] = n1 * b1[2] + n2 * b2[2] + n3 * b3[2]
        cart[3,i] = n1 * b1[3] + n2 * b2[3] + n3 * b3[3]
    end
    cart        
end
function get_n_plane_waves(self::GVector, ecut::Float64, k::Vector{Float64}) ::Int32
    nfft = size(self.cart)[2]
    npw = 0
    for i in 1:nfft
        if norm(self.cart[:,i] + k) <= sqrt(2ecut)
            npw += 1
        end
    end
    npw    
end
function get_plane_waves(self::GVector, ecut::Float64, k::Vector{Float64},
                         npw::Int32) ::Vector{Int32}
    nfft = size(self.cart)[2]
    gs = zeros(npw)
    ipw = 0
    for i in 1:nfft
        if norm(self.cart[:,i] + k) <= sqrt(2ecut)
            ipw += 1
            gs[ipw] = i
        end
    end
    gs
end
