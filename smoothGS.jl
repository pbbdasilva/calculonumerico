using LinearAlgebra
include("createb.jl")
function smoothGS(n)
    w = calculatew(n)
    d = repeat([2], n)
    dl = repeat([-1], n-1)
    A = Tridiagonal(dl, d, dl)
    b = zeros(n)
    b[1] = 1
    # solucao exata
    x = [1-k/(n+1) for k = 1:n]
    Xk0 = zeros(n)
    Xk1 = zeros(n)
    Dk1 = 1
    L = LowerTriangular(A)
    U = A - L
    iters = 0
    TOL = 10^(-8)
    while Dk1 > TOL
        if iters > 9999
            println("fracassou")
            return nothing
        end
        for k=1:n
            Xk0[k] = Xk1[k]
            Xk1[k] = (1-w)*Xk0[k] + w*(b[k] - dot(A[k,:],Xk1) + A[k,k]*Xk1[k])/A[k,k]
        end
        iters += 1
        Dk1 = maximum(abs.(Xk1 - Xk0))
    end
    println(iters)
    Erro = 100*norm(Xk1 - x)/norm(x)
    return Erro
end
