%Author: Dominic (Zhongda) Huang
%Date: 2021.08.24
%Input: evaluated point list Point_List, current best solution x_best_index,
%search radius rou
%Output: the list of evaluated point-value pair taht is within the ball of
%radius rou at x_best.


function M = pointsWithinBall(Point_List, x_best_index, rou, n)

    %Initialization
    M(1:n) = Evaluated_Point;  %initilize the list of point-value pair
    k = 0; %qualified points counter
    x_best = Point_List(x_best_index).Point;
    
    for i = 1:n
        x_current = Point_List(i).Point;
        f_current = Point_List(i).Value;
        if (norm(x_current - x_best) <= rou && ~isinf(f_current) && ~isnan(f_current))
            k = k + 1;
            M(k) = Point_List(i);
        end   
    end       
    
    M = M(1:k);
end