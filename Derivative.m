function jacobi = Derivative(op, sp, aps, eps, dt, index)
%Derivative: Numerical derivative by calculateing (f(x+dt)-f(x-dt))/2dt 

    % input: op (1 x 3)
    %        sp (1 x 3)

ap_count = length(aps);

if (index < ap_count+1)  % aps
   cur_aps = aps;
   cur_aps(index)=aps(index)+dt;
   f2 = ObsFunction(op, sp, cur_aps, eps);
   
   cur_aps(index)=aps(index)-dt;
   f1 = ObsFunction(op, sp, cur_aps, eps);
else   % eps
   cur_eps = eps;
   eps_index = mod(index-ap_count-1,6)+1;
   cur_eps(eps_index)=eps(eps_index)+dt;
   f2 = ObsFunction(op, sp, aps, cur_eps); % 3 x 1
   
   cur_eps(eps_index)=eps(eps_index)-dt;
   f1 = ObsFunction(op, sp, aps, cur_eps); % 3 x 1
end

jacobi = (f2-f1)/(2*dt); % 3 x 1

end

