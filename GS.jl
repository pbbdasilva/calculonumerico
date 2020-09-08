using LinearAlgebra
function solvegs(n, w)
    d = repeat([2], n)
    dl = repeat([-1], n-1)
    A = Tridiagonal(dl, d, dl)
    b = zeros(n)
    b[1] = 1
    # solucao exata
    x = [1-k/(n+1) for k = 1:n]
    Xk0 = zeros(1, n)
    Xk1 = zeros(1, n)
    Dk1 = max(abs(Xk1 - Xk0))
    iters = 0
    TOL = 10^(-8)
    while Dk1 > TOL
        if iter > 9999
            println("fracassou")
            return nothing
        end
        for i = 1:n
            Xk0[i] = Xk1[i]
            Xk1[i] = (b[i] - dot(A[i,:],Xk0) + A[i,i])/A[i,i]
        end
        K += 1
        Dk1 = max(abs(Xk1 - Xk0))
    end
    Erro = 100*norm(Xk1 - x)/norm(x)
    return Erro
end
