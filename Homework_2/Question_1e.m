%{
This script generates error vs work data for all methods and plots it.
%}
% Import ODE Solvers
addpath('ODE_Solvers/')
options = optimset('Display','off');

% Set global variable to count the number of function evaluations
global f_eval_count;
f_eval_count = 0;

% ODE parameters
f = @example_func_with_counter;
t0 = 0;
tf = 1;
y0 = 0;

% Exact solution
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_exact, y_exact] = ode45(f, [t0, tf], y0,opts);
exact_solution = y_exact(end);

% Stepsizes to test
h_values = 1./(100:50:500);

% List of solvers
solver_names = {'explicitEuler', 'implicitEuler', 'implicitMidpoint', 'trapezoidal', 'adamsBashford2', 'rungekutta2', 'rungekutta4'};

% Measure and plot error vs work relation
figure;
axes('XScale', 'log', 'YScale', 'log')

hold on;
for solver_idx = 1:length(solver_names)
    % compute error vs work for each solver
    solver_name = solver_names{solver_idx};

    % arrays to store error vs work data
    error_values = zeros(size(h_values));
    work_values = zeros(size(h_values));
    
    for idx = 1:length(h_values)
        % Compute error vs work for each time step choice
        h = h_values(idx);
        
        % Rest counter
        f_eval_count = 0;

        % Call the respective solver function
        [t, y] = feval(solver_name, f, t0, tf, y0, h); 

        % Calculate
        numeric_solution = y(end);
        global_error = abs(exact_solution - numeric_solution);

        % Store error and work
        error_values(idx) = global_error;
        work_values(idx) = f_eval_count;
    end
    % plot error vs work
    hold on;
    loglog(error_values,work_values,'-','DisplayName',solver_name)

end
legend('show')
title('Error vs Work', 'FontSize', 12);
xlabel('e');
ylabel('M');
saveas(gcf, 'Question_1_plts/Error_vs_Work.png');

