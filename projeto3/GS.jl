using LinearAlgebra
function solvegs(n)
    d = repeat([2], n)
    dl = repeat([-1], n-1)
    A = Tridiagonal(dl, d, dl)
    b = zeros(n)
    b[1] = 1
    # solucao exata
    x = [1-k/(n+1) for k = 1:n]
    x = reshape(x, (1,n))
    Xk0 = zeros(1, n)
    Xk1 = zeros(1, n)
    Dk1 = 1
    iters = 0
    TOL = 10^(-8)
    while Dk1 > TOL
        if iters > 9999
            println("fracassou")
            return nothing
        end
        for i = 1:n
            Xk0[i] = Xk1[i]
            Xk1[i] = (b[i] - dot(A[i,:],Xk0) + A[i,i]*Xk1[i])/A[i,i]
        end
        iters += 1
        Dk1 = maximum(abs.(Xk1 - Xk0))
        println("Após "*string(iters)*" iteraçoes:")
    end
    Erro = 100*norm(Xk1 - x)/norm(x)
    return Erro, iters
end
