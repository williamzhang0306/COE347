%{
This script plot the solutions from all methods for a single ODE on the
same plot. This is intended to demonstrate that each method works properly
%}

% Import ODE Solvers
addpath('ODE_Solvers/')
options = optimset('Display','off');


demo_f = @(t,y) -50 * (y-cos(t));

t0 = 0;
tf = 1;
y0 = 0;
h = .01;

[t_xact, y_xact] = ode45(demo_f, [0,1], 0);
[t, y] = explicitEuler(demo_f, t0, tf, y0, h);
[t2, y2] = implicitEuler(demo_f, t0, tf, y0, h);
[t3, y3] = implicitMidpoint(demo_f, t0, tf, y0, h);
[t4, y4] = trapezoidal(demo_f, t0, tf, y0, h);
[t5, y5] = adamsBashford2(demo_f, t0, tf, y0, h);
[t6, y6] = rungekutta2(demo_f, t0, tf, y0, h);
[t7, y7] = rungekutta4(demo_f, t0, tf, y0, h);

%plot
figure;
hold on
plot(t_xact, y_xact, '-', 'DisplayName', 'Exact');
plot(t, y, '-.','DisplayName', 'Explicit Euler');
plot(t2, y2, '-.','DisplayName', 'Implicit Euler');
plot(t3, y3, '-.','DisplayName', 'Implicit Midpoint');
plot(t4, y4, '-.','DisplayName', 'Trapezoid');
plot(t5, y5, '-.','DisplayName', 'Adams Bashforth 2');
plot(t6, y6, '-.','DisplayName', 'Runge Kutta 2');
plot(t7, y7, '-.','DisplayName', 'Runge Kutta 4');
legend('show', 'location', 'southeast');

xlabel('x');
ylabel('y(x)');
title('Various Solutions to y''(x) = -50(y - cos(x)) with y(0) = 0, h = 0.01', 'FontSize',12);
grid on;

% save
f = gcf;
exportgraphics(f,'Question_1_plots/All_Solvers.png','Resolution',300)