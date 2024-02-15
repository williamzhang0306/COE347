function [t, y] = rungekutta2(f, t0, tf, y0, h)
% Runge Kutta 2 ODE solver implemention
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
        t2 = t(i) + h/2;
        y2 = y(i) + (h/2)*f(t(i),y(i));
        y(i + 1) = y(i) + h * f(t2, y2);
    end
end