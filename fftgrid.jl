using LinearAlgebra

struct FFTGrid
    n1 :: Int32
    n2 :: Int32
    n3 :: Int32
end
function FFTGrid(latt::Lattice, ecutrho::Float64)
    gmax = sqrt(2ecutrho)
    n1 = ceil((2gmax * norm(latt.v1) / 2π))
    n2 = ceil((2gmax * norm(latt.v2) / 2π))
    n3 = ceil((2gmax * norm(latt.v3) / 2π))
    FFTGrid(n1, n2, n3)
end
function Base.show(io::IO, fftgrid::FFTGrid)
    print(io, "FFTGrid: ", fftgrid.n1, " ", fftgrid.n2, " ", fftgrid.n3)
end
function fft_left_end(n::Int32) ::Int32
    if n % 2 == 0
        -(n - 2) / 2
    else
        -(n - 1) / 2
    end
end
function fft_right_end(n::Int32) ::Int32
    if n % 2 == 0
        n / 2
    else
        (n - 1) / 2
    end
end

