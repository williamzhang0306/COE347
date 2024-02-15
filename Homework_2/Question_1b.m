%{
This script creates a seperate plot for each method and shows that it
converges towards the exact solution
%}

% Import ODE Solvers
addpath('ODE_Solvers/')

% ODE and parameters
demo_f = @(t, y) -50 * (y - cos(t));
t0 = 0;
tf = 1;
y0 = 0;

% Define step sizes
h_values = [0.02, 0.01, 0.005, 0.001];

% Loop over solvers
solver_names = {'explicitEuler', 'implicitEuler', 'implicitMidpoint', 'trapezoidal', 'adamsBashford2', 'rungekutta2', 'rungekutta4'};

for solver_idx = 1:length(solver_names)
    solver_name = solver_names{solver_idx};

    % Create a new figure for each solver
    figure;
    hold on;

    % Loop over step sizes
    for h_idx = 1:length(h_values)
        h = h_values(h_idx);

        % Call the respective solver function
        [t, y] = feval(solver_name, demo_f, t0, tf, y0, h);

        % Plot the solution for the current step size
        plot(t, y, '-.', 'DisplayName', ['h = ', num2str(h)]);
    end

    % Call ode45 for exact solution
    [t_exact, y_exact] = ode45(demo_f, [t0, tf], y0);
    
    % Plot the 'exact' solution
    plot(t_exact, y_exact, '-', 'DisplayName', 'Exact');

    % Add labels and legend
    xlabel('x');
    ylabel('y(x)');
    title(['Solution for ', solver_name, ' with Different Step Sizes'], 'FontSize', 12);
    legend('show', 'location', 'southeast');

    % Save the plot
    saveas(gcf, ['Question_1_plots',solver_name, '_converging.png']);
    close(gcf); % Close the figure to avoid accumulating plots
end