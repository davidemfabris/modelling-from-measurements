function rhs=lorenz_rhs(t,x,dummy,sigma,beta,rho)

rhs=[sigma*(-x(1)+x(2))
     -x(1)*x(3)+rho*x(1)-x(2)
     x(1)*x(2)-beta*x(3)];