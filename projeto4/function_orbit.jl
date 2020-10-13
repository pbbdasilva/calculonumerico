using Polynomials
function orbit()
    r0 = 6870;
    r1 = 6728;
    r2 = 6675;
    theta0 = -pi/6;
    theta1 = 0;
    theta2 = pi/6;
    y = [r0,r1,r2];
    x = [theta0,theta1,theta2];
    poly = fit(x, y);
    return poly
end
