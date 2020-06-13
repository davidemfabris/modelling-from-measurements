function rhs=lotka_volterra_rhs(t,x,dummy,Mat)

rhs=Mat*[x(1);x(2);x(1).*x(2)];