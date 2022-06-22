% Author: Dominic Huang
% Date: 2021.09.08
% Linear Step for rqlif method
% Input:
%       - Objective funcction f
%       - Current best point x_best
%       - Step Length
% Output: Evaluated point list

function x_L = linearStep03(f, Point_List, x_best_index)

x_left = Point_List(x_best_index(end - 1));

x_mid = Point_List(x_best_index(end));

x_right = 2 * x_mid.Point - x_left.Point;

x_L(1) = Evaluated_Point;
x_L(1).Point = x_right;
x_L(1) = eval(f, x_L(1));
end