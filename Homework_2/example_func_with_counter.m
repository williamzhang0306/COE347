function f = example_func_with_counter(t,y)
% returns -50 * (y-cos(t)) )The example ODE)
% Also increases a global counter for the number of times f(t,y) is
% evaluated.

    f = -50 * (y-cos(t));

    % increase global variable counter
    global f_eval_count;
    f_eval_count = f_eval_count + 1;

end

