function [t, y] = rungekutta4(f, t0, tf, y0, h)
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
        k1 = f(t(i), y(i));
        k2 = f(t(i)+h/2, y(i)+(h/2)*k1);
        k3 = f(t(i)+h/2, y(i)+(h/2)*k2);
        k4 = f(t(i)+h, y(i)+h*k3);
        y(i + 1) = y(i) + h*(k1/6 + k2/3 + k3/3 +k4/6);
    end
end