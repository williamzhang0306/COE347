%{
This script calculates the global error for each method at a range of step
sizes and produces a loglog plot.
%}
% Import ODE Solvers
addpath('ODE_Solvers/')

% Define ODE and parameters
demo_f = @(t, y) -50 * (y - cos(t));
t0 = 0;
tf = 1;
y0 = 0;

% Get exact solution
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_exact, y_exact] = ode45(demo_f, [t0, tf], y0,opts);
exact_solution = y_exact(end);

% Stepsizes to test
%h_values = [0.02, 0.01, 0.005, 0.001];

one_over_h = 100:50:500;
h_values = 1./one_over_h;

% List of solvers
solver_names = {'explicitEuler', 'implicitEuler', 'implicitMidpoint', 'trapezoidal', 'adamsBashford2', 'rungekutta2', 'rungekutta4'};

figure;
axes('XScale', 'log', 'YScale', 'log')
hold on;
for solver_idx = 1:length(solver_names)

    solver_name = solver_names{solver_idx};

    % Create a new figure for each solver


    error_values = zeros(size(h_values));

    for idx = 1:length(h_values)
        h = h_values(idx);

        % Call the respective solver function
        [t, y] = feval(solver_name, demo_f, t0, tf, y0, h); 

        % Calculate and save error
        numeric_solution = y(end);
        global_error = abs(exact_solution - numeric_solution);
        error_values(idx) = global_error;

    end
    plot(one_over_h,error_values,'o-','DisplayName',solver_name)
end

legend('show','Location','southwest');
title('Global Error vs Step Size for Various Numerical Methods', 'FontSize', 12);
xlabel('1/h');
ylabel('Error')
saveas(gcf, 'Question_1_plots/Error_Plot.png');
