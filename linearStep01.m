% Author: Dominic Huang
% Date: 2021.09.08
% Linear Step for rqlif method
% Input: 
%       - Objective funcction f
%       - Current best point x_best
%       - Step Length
% Output: Evaluated point list

function x_L = linearStep01(f, x_best, step)
    x_L(1) = Evaluated_Point;
    x_L(1).Point = x_best.Point - x_best.Gradient / norm(x_best.Gradient) * step;
    x_L(1) = eval(f, x_L(1));
end