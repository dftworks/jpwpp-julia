function simpson_log(y :: Vector{Float64}, r :: Vector{Float64}) :: Float64
    n = size(y)[1]
    dx = log((r[n]/r[1])) / (n-1.0)
    s = 0.0
    for i in 1:2:n-1
        s += 1.0 * y[i] * r[i]
        s += 4.0 * y[i+1] * r[i+1]
        s += 1.0 * y[i+2] * r[i+2]
    end
    s * dx / 3.0
end
