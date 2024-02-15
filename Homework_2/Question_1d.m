%{
This script fits a power law to error vs stepsize data and also produces
plots for visualizing
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

h_values = 1./(100:20:500);
one_over_h = 1./h_values;

% List of solvers
solver_names = {'explicitEuler', 'implicitEuler', 'implicitMidpoint', 'trapezoidal', 'adamsBashford2', 'rungekutta2', 'rungekutta4'};

% Power law to fite e(h) to
powerLawModel = fittype('C * x^a', 'coefficients', {'C', 'a'}, 'independent', 'x');


for solver_idx = 1:length(solver_names)
    solver_name = solver_names{solver_idx};

    % Create a new figure for each solver
    figure;
    hold on;

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
    
    % Log-transform data for linear fit
    log_h_values = log(h_values);
    log_error_values = log(error_values);

    % Perform linear fit in log space
    fitResult = polyfit(log_h_values, log_error_values, 1);
    alpha = fitResult(1);
    C = exp(fitResult(2));
    disp([solver_name, ' alpha = ', num2str(alpha)])

    % Plot original data and linear fit
    plot(h_values, error_values, 'o', 'DisplayName', [solver_name,' global error']);
    plot(h_values, exp(polyval(fitResult, log_h_values)), 'DisplayName', [num2str(C),'*x ^{',num2str(alpha),'}']);
    

    % plot options
    title(['Error Order Analysis: ', solver_name], 'FontSize', 20);
    xlabel('h');
    ylabel('error')
    xlim([0 1.2*max(h_values)])
    ylim([0, 1.2*max(error_values)])
    legend('show','FontSize', 16, 'Location','northwest')
    
    % save ooptions
    saveas(gcf, ['Question_1_plots/',solver_name, '_ErrorOrder.png']);
    close(gcf); % Close the figure to avoid accumulating plots
end
