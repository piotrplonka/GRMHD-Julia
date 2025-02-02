function lu_decomposition(A::Array{Float64,2})
    n = size(A, 1)
    L = zeros(Float64, n, n)
    U = copy(A)
    p = collect(1:n)
    for k in 1:n
        pivot = k
        maxval = abs(U[k, k])
        for i in k+1:n
            if abs(U[i, k]) > maxval
                maxval = abs(U[i, k])
                pivot = i
            end
        end
        if abs(U[pivot, k]) < 1e-12
            error("Macierz jest osobliwa lub źle uwarunkowana.")
        end

        if pivot != k
            U[k, :], U[pivot, :] = U[pivot, :], U[k, :]
            if k > 1
                L[k, 1:k-1], L[pivot, 1:k-1] = L[pivot, 1:k-1], L[k, 1:k-1]
            end
            p[k], p[pivot] = p[pivot], p[k]
        end

        for i in k+1:n
            L[i, k] = U[i, k] / U[k, k]
            for j in k:n
                U[i, j] -= L[i, k] * U[k, j]
            end
        end

        L[k, k] = 1.0
    end

    return L, U, p
end

function forward_substitution(L::Array{Float64,2}, b::Vector{Float64})
    n = length(b)
    y = zeros(Float64, n)
    for i in 1:n
        sum = 0.0
        for j in 1:i-1
            sum += L[i, j] * y[j]
        end
        y[i] = b[i] - sum
    end
    return y
end

function backward_substitution(U::Array{Float64,2}, y::Vector{Float64})
    n = length(y)
    x = zeros(Float64, n)
    for i in n:-1:1
        sum = 0.0
        for j in i+1:n
            sum += U[i, j] * x[j]
        end
        x[i] = (y[i] - sum) / U[i, i]
    end
    return x
end

function invert_matrix_lu_manual(A::Array{Float64,2})
    n = size(A, 1)
    L, U, p = lu_decomposition(copy(A))
    invA = zeros(Float64, n, n)
    for i in 1:n
        e = zeros(Float64, n)
        e[i] = 1.0
        b = zeros(Float64, n)
        for j in 1:n
            b[j] = e[p[j]]
        end
        y = forward_substitution(L, b)
        x = backward_substitution(U, y)
        invA[:, i] = x
    end
    return invA
end


function solve_lu(A::Array{Float64,2}, b::Vector{Float64})
    L, U, p = lu_decomposition(copy(A))
    n = length(b)
    bp = zeros(Float64, n)
    for i in 1:n
        bp[i] = b[p[i]]
    end
    y = forward_substitution(L, bp)
    x = backward_substitution(U, y)
    return x
end
