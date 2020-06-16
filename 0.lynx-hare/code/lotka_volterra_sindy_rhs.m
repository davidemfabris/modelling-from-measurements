function rhs=lotka_volterra_sindy_rhs(t,x,dummy,A_handle, param)

rhs =  param*A_handle(x(1), x(2))';
     