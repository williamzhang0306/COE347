function [t, y] = adamsBashford2(f, t0, tf, y0, h)
% Explicity Euler ODE solver implemention
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

    % use rk2 to start
    [sub_t, sub_y] = rungekutta2(f,t0,t0+h,y0,h);
    t(2) = sub_t(2);
    y(2) = sub_y(2);
    
    for i = 2:(length(t) - 1)
        y(i + 1) = y(i) + h * (1.5*f(t(i),y(i)) - .5*f(t(i-1),y(i-1)));
    end
end