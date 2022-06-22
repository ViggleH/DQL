# DQL
The MATLAB implementation of the DQL method and Smart DQL method. The DQL method is a hybrid derivative free optimization method. It has mathematically proven local convergence. The methods are part of my thesis work as follows:
Huang, Z. (2022). A hybrid direct search and model-based derivative-free optimization method with dynamic decision processing.

[Point_List, x_best_index, stop, time, step, n, k] = DQL(f, x_0, step, eps_max, eps_min, eps_grad, max_calls, direct_switch, quad_switch, linear_switch, dispOption)
[Point_List, x_best_index, stop, time, step, n, k] = sDQL(f, x_0, step, eps_max, eps_min, eps_grad, max_calls, dispOption)


Input:
   - Objective Function: f
   - Search Start Point: x_0
   - Search Step Length: step
   - Disired Accuracy: eps_max
   - Numerical Error Safeguard: eps_min
   - Disired gradient: eps_grad
   - Maximum Function Calls: max_calls
   - Maximum search stage: max_search
   - Direct Swtich:
       - 0: Disable direct step
       - 1: Random rotation
       - 2: Smart direct step (should be use with quadratic switch 3 and linear switch 4)
   - Quadratic Swtich:
       - 0: Disable quadratic step
       - 1: Quadratic step with quadratic model
       - 2: Quadratic step with Newton Step & trust region
       - 3: Smart quadratic step
   - Linear Switch:
       - 0: Disable linear step
       - 1: Search on the steepest descent direction with fixed step
       length
       - 2: The algorithm tries to construct a bracket on the steepest
       descent direction and searches with safeguarded sequential
       quadratic fitting algorithm at 1 and 1/2 step.
       - 3: Search on descent direction from last x_best
       - 4: Smart Linear Step
    -dispOption:
      -0: disable display
      -1: display the estimated remaining time.
\{Output}:
   - Evaluated Point_List
   - The indexes of best solutions.
   - Stop Indicator: - 1: optimum solution found
                     - 0: optimum solution not found
                     - -1: optimum solution not found because the search step does not meet condition.
                     - -2: optimum solution not found because the search step and gradient does not meet condition.

