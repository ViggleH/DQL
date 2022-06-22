% Author: Dominic Huang
% Date: 2021.09.08
% Linear Step for rqlif method
% Input:
%       - Objective funcction f
%       - Current best point x_best
%       - Step Length
% Output: Evaluated point list

function x_Linear = linearStep05(f, Point_List, x_best_index, step, eps1)
%Compare best point distance to step
descent_length = 0.5 * norm(Point_List(x_best_index(end)).Point - Point_List(x_best_index(end - 1)).Point);
n = 3; %number of iternation


if step > eps1
    if descent_length > step
        
        x_left = Point_List(x_best_index(end - 1));
        
        x_mid = Point_List(x_best_index(end));
        
        x_right = 2 * x_mid.Point - x_left.Point;
        
        x_Linear(1) = Evaluated_Point;
        x_Linear(1).Point = x_right;
        x_Linear(1) = eval(f, x_Linear(1));
    else
        x_Linear(1) = Evaluated_Point;
        x_Linear(1).Point = Point_List(x_best_index(end)).Point - Point_List(x_best_index(end)).Gradient / norm(Point_List(x_best_index(end)).Gradient) * step;
        x_Linear(1) = eval(f, x_Linear(1));
    end
else
        x_best = Point_List(x_best_index(end));
        x_right = Evaluated_Point;
        x_right.Point = x_best.Point - x_best.Gradient / norm(x_best.Gradient) * step;
        x_right = eval(f, x_right);
        
        x_mid = Evaluated_Point;
        x_mid.Point = x_best.Point - x_best.Gradient / norm(x_best.Gradient) * step * 1/2;
        x_mid = eval(f, x_right);
        x_Linear = [x_mid, x_right];
        if(x_mid.Value < min([x_best.Value, x_right.Value]))
            %A bracket is obtained, Line Search
            
            %Parametrise f
            F = @(a) f(x_best.Point - a * (step / norm(x_best.Gradient)) * x_best.Gradient );
            
            a_l = 0;
            a_0 = 1/2;
            a_r = 1;
            F_l = x_best.Value;
            F_0 = x_mid.Value;
            F_r = x_right.Value;
            
            [~, ~, XF_Line] = bracketLineSearch(F, a_l, F_l, a_0, F_0, a_r, F_r, n);
            for i = 1:size(XF_Line, 2)
                x_Linear(i+2).Point = x_best.Point + XF_Line(1,i) * (step / norm(x_best.Gradient)) * x_best.Gradient;
                x_Linear(i+2).Value = XF_Line(2,i);
                x_Linear(i+2).Type = 0;
            end
        end
end
end