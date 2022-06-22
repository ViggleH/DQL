%Author: Dominic (Zhongda) Huang
%Date: 2022.05.02
%(Smart) DQL method for a box-constrained black-box function
%Input:
%   - Objective Function: f
%   - Search Start Point: x_0
%   - Search Step Length: step
%   - Disired Accuracy: eps_max
%   - Numerical Error Safeguard: eps_min
%   - Disired gradient: eps_grad
%   - Maximum Function Calls: max_calls
%   - Maximum search stage: max_search
%   - Direct Swtich:
%       - 0: Disable direct step
%       - 2: S should be use with quadratic switch 3 and linear switch 4
%   - Quadratic Swtich:
%       - 0: Disable quadratic step
%       - 1: Quadratic step with quadratic model
%       - 2: Quadratic step with Newton Step & trust region
%       - 3: S
%   - Linear Switch:
%       - 0: Disable linear step
%       - 1: Search on the steepest descent direction with fixed step
%       length
%       - 2: The algorithm tries to construct a bracket on the steepest
%       descent direction and searches with safeguarded sequential
%       quadratic fitting algorithm at 1 and 1/2 step.
%       - 3: Search on descent direction from last x_best
%       - 4: S
%Output:
%   - Evaluated Point_List
%   - The indexes of best solutions.
%   - Stop Indicator: - 1: optimum solution found
%                     - 0: optimum solution not found

function [Point_List, x_best_index, stop, time, step, n, k] = DQL(f, x_0, step, eps_max, eps_min, eps_grad, max_calls, direct_switch, quad_switch, linear_switch, dispOption)

%Initialization
tic;
dim = size(x_0, 1); %Dimension of the fucntion input
Point_List(1:max_calls) = Evaluated_Point;
q = 0.3;    %Update parameter
rou = 3 * step; %Search radius in quadratic step
t_min = 2 * dim; % minimum points required to initiate quadratic step
max_search = max_calls * 0.9;


Point_List(1).Point = x_0;
Point_List(1) = eval(f, Point_List(1)); %Evaluate f(x_0);
x_best_index = zeros(1, max_calls);
x_best_index(1) = 1;
b = 1;

n = 1;  %Evaluation counter. 1 for f(x_0).
k = 1;  %Iteration counter
D1 = [eye(dim), -eye(dim)];  %Default search directions
stop = 0;   %Stop Indicator

tau = 0; %Regularization Parameter;

FlagDirect = 0;
FlagQuadratic = 0;
FlagLinear = 0;
x_L = [];
x_center_index = zeros(1, max_calls);

while(stop == 0)
    if nargin > 10
        if(dispOption == 1)
            disp('Optimizing Progress: ' + string(n/max_calls * 100) + '%; Elapsed Time:' + string(toc) + 's; ' + 'Estimate Remaining: ' + string((max_calls-n) * toc/n) + 's.');
        end
    end
    %Update the indexes for initial x_best
    x_center_index(k) = x_best_index(b);
    
    %Direct Step
    if(direct_switch ~= 0)
        
        %Determine the rotation matrix rot
        switch(direct_switch)
            case 1
                if(isEven(k))
                    [ra, rb] = randomHyperplane(dim);
                    rot = rotMatrix(pi * 2 * rand(1), ra, rb);
                else
                    rot = eye(dim);
                end
                %{
            case 2
                if(Point_List(x_best_index(b)).Type >= 1)
                    rot = alignMatrix(D1, -Point_List(x_best_index(b)).Gradient);
                else
                    rot = eye(dim);
                end
            case 3
                if(Point_List(x_best_index(b)).Type >= 1)
                    d1 = ones(dim, 1);
                    rot = alignMatrix(d1, -Point_List(x_best_index(b)).Gradient);
                else
                    rot = eye(dim);
                end
            case 4
                if FlagQuadratic == 1
                        if FlagNewton == 0 && ~isempty(quad_grad)
                            descent = -quad_grad;
                            if norm(descent, inf) >= eps3
                                rot = alignMatrix(D1, descent);
                            else
                                rot = eye(dim);
                            end
                        else
                            rot = eye(dim);
                        end
                else
                        rot = eye(dim);
                end
                %}
            case 2
                if FlagDirect == 0
                    if FlagQuadratic == 1
                        if FlagNewton == 0 && ~isempty(quad_grad)
                            descent = -quad_grad;
                            if norm(descent, inf) >= eps_grad
                                rot = alignMatrix(D1, descent);
                            else
                                rot = eye(dim);
                            end
                        else
                            rot = eye(dim);
                        end
                    else
                        if FlagLinear == 1
                            descent = Point_List(x_best_index(b)).Point - Point_List(x_center_index(k-1)).Point;
                            rot = alignMatrix(D1, descent);
                        else
                            rot = eye(dim);
                            for i = 1 : size(x_L, 2)
                                if(sum(isnan(x_L(i).Point)) == 0 && ~isempty(x_L(i).Point))
                                    d1 = ones(dim, 1);
                                    ascent = x_L(i).Point - Point_List(x_center_index(k)).Point;
                                    rot = alignMatrix(d1, ascent);
                                    break;
                                end
                            end
                        end
                    end
                end
        end
        
        %Initialize search direction D1
        D = rot * D1;
        
        
        %Initialize the indicators for the results from each step: 0 means fail; 1 means success.
        FlagDirect = 0;
        FlagQuadratic = 0;
        FlagLinear = 0;
        FlagNewton = 0;
        
        %Search the solution on directions defined by D at the fixed step length
        Direct_Step_List = directStep(f, Point_List(x_best_index(b)), step, D);
        
        %Update center gradient of x_best and the points in direct step
        Point_List(x_best_index(b)) = grad(Point_List(x_best_index(b)), Direct_Step_List);
        
        %Append the direct step list to evaluated list, update x_best if there
        %is a better point.
        for i = 1 : size(Direct_Step_List, 2)
            Point_List(n + i) = Direct_Step_List(i);
            if(k <= max_search)
                if Point_List(n + i).Value < Point_List(x_best_index(b)).Value
                    b = b + 1;
                    x_best_index(b) = n + i;
                    FlagDirect = 1;
                end
            else
                if Point_List(n + i).Value < Point_List(x_best_index(b)).Value - step^2
                    b = b + 1;
                    x_best_index(b) = n + i;
                    FlagDirect = 1;
                end
            end
        end
        n = n + size(Direct_Step_List, 2);
        
    end
    
    if(k > max_search)
        quad_switch = 0;
        linear_switch = 0;
    end
    
    %Check if the stopping condition is met
    if(Point_List(x_best_index(b)).Type >= 1 && norm(Point_List(x_best_index(b)).Gradient, inf) < eps_grad && step < eps_max || n + 2*dim + 6 > max_calls)
        stop = 1;
        break;
    end
    
    %Quad Step
    if(quad_switch ~= 0)
        
        switch (quad_switch)
            case 1
                %Find the evaluated points within the ball of radius rou at
                %x_best, store the list as M.
                M = pointsWithinBall(Point_List, x_best_index(b), rou, n);
                
                %Initiate quadratic step if enough points are obtained
                if(size(M,2) >= t_min)
                    x_Q = quadraticStep(f, Point_List(x_best_index(b)), M, tau, step);
                    %If an improvement is found, update x_best, f_best and
                    %the indicator.
                    if(sum(isnan(x_Q.Point), 'all') + sum(isinf(x_Q.Point), 'all') == 0  && ~isempty(x_Q.Point))
                        x_Q = eval(f, x_Q);
                        Point_List(n + 1) = x_Q;
                        if x_Q.Value < Point_List(x_best_index(b)).Value
                            b = b + 1;
                            x_best_index(b) = n + 1;
                            FlagQuadratic = 1;
                        end
                        n = n + 1;
                    end
                end
            case 2
                M = pointsWithinBall(Point_List, x_best_index(b), rou, n);
                h = 0;
                pointsWithGrad(1:size(M,2)) = Evaluated_Point;
                for i = 1:size(M,2)
                    if(M(i).Type == 1)
                        h = h + 1;
                        pointsWithGrad(h) = M(i);
                    end
                end
                pointsWithGrad = pointsWithGrad(1:h);
                
                if(Point_List(x_best_index(b)).Type == 1 && h~= 0)
                    Point_List(x_best_index(b)) = hes(Point_List(x_best_index(b)), pointsWithGrad);
                    
                    %check the definiteness of the Hessian
                    if(Point_List(x_best_index(b)).Type == 2)
                        d = eig(Point_List(x_best_index(b)).Hessian);
                        if(all(d > 0))
                            x_Q = Evaluated_Point;
                            x_Q.Point = Point_List(x_best_index(b)).Point - Point_List(x_best_index(b)).Hessian \ Point_List(x_best_index(b)).Gradient;
                            
                            if(sum(isnan(x_Q.Point), 'all') + sum(isinf(x_Q.Point), 'all')  == 0 && ~isempty(x_Q.Point))
                                x_Q = eval(f, x_Q);
                                FlagNewton =  1;
                                Point_List(n + 1) = x_Q;
                                if x_Q.Value < Point_List(x_best_index(b)).Value
                                    b = b + 1;
                                    x_best_index(b) = n + 1;
                                    FlagQuadratic = 1;
                                end
                                n = n + 1;
                            end
                        else
                            x_Q = Evaluated_Point;
                            x_temp = trust(Point_List(x_best_index(b)).Gradient, Point_List(x_best_index(b)).Hessian, step*2);
                            x_Q.Point = Point_List(x_best_index(b)).Point + x_temp;
                            
                            
                            if(sum(isnan(x_Q.Point), 'all') + sum(isinf(x_Q.Point), 'all')  == 0 && ~isempty(x_Q.Point) && isreal(x_Q.Point))
                                quad_grad = (Point_List(x_best_index(b)).Hessian + Point_List(x_best_index(b)).Hessian') * x_temp + Point_List(x_best_index(b)).Gradient;
                                x_Q = eval(f, x_Q);
                                Point_List(n + 1) = x_Q;
                                if x_Q.Value < Point_List(x_best_index(b)).Value
                                    b = b + 1;
                                    x_best_index(b) = n + 1;
                                    FlagQuadratic = 1;
                                end
                                n = n + 1;
                            end
                            
                            
                        end
                    end
                end
                
            case 3
                
                M = pointsWithinBall(Point_List, x_best_index(b), rou, n);
                
                h = 0;
                for i = 1:size(M,2)
                    if(M(i).Type == 1)
                        h = h + 1;
                        pointsWithGrad(h) = M(i);
                    end
                end
                
                if(Point_List(x_best_index(b)).Type == 1)
                    Point_List(x_best_index(b)) = hes(Point_List(x_best_index(b)), pointsWithGrad);
                end
                if(Point_List(x_best_index(b)).Type == 2)
                    %check the definiteness of the Hessian
                    d = eig(Point_List(x_best_index(b)).Hessian);
                    if(all(d > 0))
                        x_Q = Evaluated_Point;
                        x_Q.Point = Point_List(x_best_index(b)).Point - Point_List(x_best_index(b)).Hessian \ Point_List(x_best_index(b)).Gradient;
                        
                        
                        if(sum(isnan(x_Q.Point), 'all') + sum(isinf(x_Q.Point), 'all') == 0  && ~isempty(x_Q.Point))
                            x_Q = eval(f, x_Q);
                            FlagNewton =  1;
                            Point_List(n + 1) = x_Q;
                            if x_Q.Value < Point_List(x_best_index(b)).Value
                                b = b + 1;
                                x_best_index(b) = n + 1;
                                FlagQuadratic = 1;
                            end
                            n = n + 1;
                        end
                    else
                        x_Q = Evaluated_Point;
                        x_temp = trust(Point_List(x_best_index(b)).Gradient, Point_List(x_best_index(b)).Hessian, step*2);
                        x_Q.Point = Point_List(x_best_index(b)).Point + x_temp;
                        
                        
                        if(sum(isnan(x_Q.Point), 'all') + sum(isinf(x_Q.Point), 'all')  == 0 && ~isempty(x_Q.Point) && isreal(x_Q.Point))
                            quad_grad = (Point_List(x_best_index(b)).Hessian + Point_List(x_best_index(b)).Hessian') * x_temp + Point_List(x_best_index(b)).Gradient;
                            x_Q = eval(f, x_Q);
                            Point_List(n + 1) = x_Q;
                            if x_Q.Value < Point_List(x_best_index(b)).Value
                                b = b + 1;
                                x_best_index(b) = n + 1;
                                FlagQuadratic = 1;
                            end
                            n = n + 1;
                        end
                    end
                    
                    
                else
                    if(size(M,2) >= t_min)
                        [x_Q, quad_grad] = quadraticStep(f, Point_List(x_best_index(b)), M, tau, step);
                        %If an improvement i s found, update x_best, f_best and
                        %the indicator.
                        if(sum(isnan(x_Q.Point), 'all') + sum(isinf(x_Q.Point), 'all')  == 0 && ~isempty(x_Q.Point))
                            x_Q = eval(f, x_Q);
                            Point_List(n + 1) = x_Q;
                            if x_Q.Value < Point_List(x_best_index(b)).Value
                                b = b + 1;
                                x_best_index(b) = n + 1;
                                FlagQuadratic = 1;
                            end
                            n = n + 1;
                        end
                    end
                    
                    
                end
        end
    end
    
    %Linear Step
    if(linear_switch ~= 0)
        if(FlagDirect == 0 && FlagQuadratic == 0)
            FlagEnoughPoints = 1;
            switch(linear_switch)
                case 1
                    Point_List(x_best_index(b)) = grad(Point_List(x_best_index(b)), Direct_Step_List);
                    if Point_List(x_best_index(b)).Type >= 1
                        x_L = linearStep01(f, Point_List(x_best_index(b)), step);
                    else
                        FlagEnoughPoints = 0;
                    end
                case 2
                    Point_List(x_best_index(b)) = grad(Point_List(x_best_index(b)), Direct_Step_List);
                    if Point_List(x_best_index(b)).Type >= 1
                        x_L = linearStep02(f, Point_List(x_best_index(b)), step);
                    else
                        FlagEnoughPoints = 0;
                    end
                case 3
                    if b >= 2
                        x_L = linearStep03(f, Point_List, x_best_index(1:b));
                    else
                        FlagEnoughPoints = 0;
                    end
                case 4
                    if b >= 2
                        x_L = linearStep04(f, Point_List, x_best_index(1:b));
                    else
                        FlagEnoughPoints = 0;
                    end
                case 5
                    if b >= 2 && Point_List(x_best_index(b)).Type >= 1
                        x_L = linearStep05(f, Point_List, x_best_index(1:b), step, eps_max);
                    else
                        FlagEnoughPoints = 0;
                    end
            end
            
            %If an improvement is found, update x_best, f_best and
            %the indicator.
            if FlagEnoughPoints == 1
                for i = 1 : size(x_L, 2)
                    if(sum(isnan(x_L(i).Point)) == 0 && ~isempty(x_L(i).Point))
                        Point_List(n + 1) = x_L(i);
                        if x_L(i).Value < Point_List(x_best_index(b)).Value
                            b = b + 1;
                            x_best_index(b) = n + 1;
                            FlagLinear = 1;
                        end
                        n = n + 1;
                    end
                end
            end
            
            
        end
    end
    
    %Update Step
    if(FlagDirect == 1)
        step = step / q;
    else
        if(FlagDirect == 0 && FlagQuadratic ==0)
            step = step * q;
        end    
    end
    rou = 2 * step;
    
    %if step is too small, terminate the algorithm
    if(step < eps_min || step > 1/eps_min)
        break;
    end
    
    k = k + 1;
end
%Trim unused entries
Point_List = Point_List(1:n);
x_best_index = x_best_index(1:b);
if stop == 0
    if step < eps_max
        stop = -1;
    end
    if norm(Point_List(x_best_index(b)).Gradient, inf) < eps_grad
        stop = -2;
    end
end


time = toc;

end
