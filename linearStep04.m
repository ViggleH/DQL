% Author: Dominic Huang
% Date: 2021.09.08
% Linear Step for rqlif method
% Input:
%       - Objective funcction f
%       - Current best point x_best
%       - Step Length
% Output: Evaluated point list

function x_Linear = linearStep04(f, Point_List, x_best_index)

x_left = Point_List(x_best_index(end - 1));

x_mid = Point_List(x_best_index(end));

x_right = Evaluated_Point;
x_right.Point = 2 * x_mid.Point - x_left.Point;
x_right = eval(f, x_right);
x_Linear(1) = x_right;

if(x_mid.Value < min([x_left.Value, x_right.Value]))
    %A bracket is obtained, Line Search
    
    %Parametrise f
    F = @(a) f(x_mid.Point + a * (x_mid.Point - x_left.Point) );
    
    a_l = -1;
    a_0 = 0;
    a_r = 1;
    F_l = x_left.Value;
    F_0 = x_mid.Value;
    F_r = x_right.Value;
    
    %number of iternation
    n = 3;
    [~, ~, XF_Line] = bracketLineSearch(F, a_l, F_l, a_0, F_0, a_r, F_r, n);
    for i = 1:size(XF_Line, 2)
        x_Linear(i+1).Point = x_mid.Point + XF_Line(1,i) * (x_mid.Point - x_left.Point);
        x_Linear(i+1).Value = XF_Line(2,i);
        x_Linear(i+1).Type = 0;
    end
end
end