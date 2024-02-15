function [t, y] = explicitEuler(f, t0, tf, y0, h)
% Explicit Euler ODE solver implemention
%{
    f : f(t,y) = y'. The derivative of y at (t,y)
    t0 : initial time
    tf : final time
    y0 : initial condition (t0,y0)
    h : step size
%}
    t = t0:h:tf;
    y = zeros(length(t), 1);
    y(1) = y0;

    for i = 1:(length(t) - 1)
        y(i + 1) = y(i) + h * f(t(i), y(i));
    end
end