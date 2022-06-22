function [Point_List, x_best_index, stop, time, step, n, k] = sDQL(f, x_0, step, eps_max, eps_min, eps_grad, max_calls, dispOption)
    [Point_List, x_best_index, stop, time, step, n, k] = DQL(f, x_0, step, eps_max, eps_min, eps_grad, max_calls, 2, 3, 4, dispOption);
end
