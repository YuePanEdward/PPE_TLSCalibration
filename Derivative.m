function jacobi = Derivative(res,dt)
%Derivative: Numerical derivative by calculateing (f(x+dt)-f(x-dt))/2dt 

% TODO

f2=funct(x+dt);
f1=funct(x-dt);

jacobi = (f2-f1)/(2*dt);

end

