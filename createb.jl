import LinearAlgebra
function calculateG(n)
    d = repeat([2], n)
    dl = repeat([-1], n-1)
    resultados = Array{Float64}(undef, 20)
    w = collect(0.1:0.1:2)
    A = Tridiagonal(dl, d, dl)
    D = Diagonal(d)
    L = LowerTriangular(A) - D
    U = UpperTriangular(A) - D
    i = 1
    for w_k in w
        M = (1/w_k)*D + L
        N = ((1-w_k)/(w_k))*D - U
        M = Matrix(M)
        N = Matrix(N)
        G = M \ N
        F = svd(G)
        resultados[i] = F.S[1]
        i += 1
    end
    index = argmin(resultados)
    return w[index], G
end
