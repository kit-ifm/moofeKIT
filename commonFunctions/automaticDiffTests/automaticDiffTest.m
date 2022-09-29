x0 = dlarray([-1,2]);
[fval,gradval] = dlfeval(@rosenbrock,x0,12)


function [y,dydx] = rosenbrock(x,temp)

y = 100*(x(2) - x(1).^2).^2 + (1 - x(1)).^2;
dydx = dlgradient(y,x);
temp2 = temp;

end