function [t, y] = implicitEuler(f, t0, tf, y0, h)
% Implicit Euler ODE solver implemention
%{
    f : f(t,y) = y'. The derivative of y at (t,y)
    t0 : initial time
    tf : final time
    y0 : initial condition (t0,y0)
    h : step size
%}
    t = t0:h:tf;
    y = zeros(length(t), length(y0));
    y(1) = y0;

    for i = 1:(length(t) - 1)
        % use a solver to find the next step
        options = optimset('Display','off'); % silence output of fsolve
        tNext = t(i+1);
        implicitEq = @(yNext) yNext - y(i) - h*f(tNext, yNext);
        y(i + 1) = fsolve(implicitEq, y(i), options);
    end
end