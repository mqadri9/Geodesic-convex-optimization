function [t,y] = forwardEulerSys(fun,tspan,y0,h)
% [t,y] = forwardEuler( fun, tspan, y0, h ) solves ODE using the forward 
% Euler method.
% tspan = [T0, TF] integrates the ODE from time T0 to time TF, with initial
% condition y0, using the forward Euler method on an equispaced grid with
% step size h. Function fun(T, Y) corresponds to f(t, y)


y = y0;
%N = (tspan(2)-tspan(1))/h;
if tspan(2) >= tspan(1)
    t = tspan(1):h:tspan(2);
else
    t = tspan(1):-h:tspan(2);
    h = -h;
end

for i = 1:length(t)-1
    y(:,i+1) = y(:,i) + h*fun(t(i),y(:,i));
    %t=[t;t(length(t))+h];
end

if tspan(2) < tspan(1)
    t = t(end:-1:1);
    y = y(end:-1:1,:);
end

end

