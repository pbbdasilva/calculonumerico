using Polynomials
function raio_min(poly)
    dp = derivative(poly);
    raiz = roots(dp);
    println(raiz);
    r_min = poly(raiz[1]);
    return r_min
end
