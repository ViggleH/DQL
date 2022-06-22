%Author: Dominic (Zhongda) Huang
%Date: 2021.08.24
%Input: objective function f, current best solution x_best,
%step length step and search direction D.
%   - D is an n by m matrix witheach colomn being the vector presentation of the search direction.
%Output: evaluated point by direct step

function Direct_Step_List = directStep(f, x_best, step, D)

    %Initialization 
    Direct_Step_List(1:size(D, 2)) = Evaluated_Point; %List of point-value pair
    best_point = x_best.Point;
    
    F(size(D, 2)).f = f;
    for i = 1:size(D, 2) - 1
        F(i).f = f;
    end
    
    parfor i = 1:size(D, 2)
        x_current = Evaluated_Point;
        x_current.Point = best_point + D(:,i)*step; %Search the point on one of the search direction with fixed step
        x_current.Value = F(i).f(x_current.Point);
        x_current.Type = 0;
        Direct_Step_List(i) = x_current; %store the point-value pair in the list
    end
    
end