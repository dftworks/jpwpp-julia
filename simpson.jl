function simpson(f::Vector{Float64}, x::Vector{Float64}) ::Float64
    n = size(f)[1]
    dx = (x[n]-x[1]) / (n-1.0)
    s = 0.0
    for i in 1:2:n-1
        s += 1.0 * f[i]
        s += 4.0 * f[i+1] 
        s += 1.0 * f[i+2] 
    end
    s * dx / 3.0
end
