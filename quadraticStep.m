%Author: Dominic (Zhongda) Huang
%Date: 2021.08.24
%Input: objective function f, current best solution x_best, constraint
%const, list of evaluated point-value pairs M, and regularization parameter
%tau.
%Output: best solution by quadratic step, x_Q.

function [x_Q, grad] = quadraticStep(f, x_best, M, tau, step)

%Initialization
grad = [];
dim = size(x_best.Point, 1);  %dimension of the problem
x_Q = Evaluated_Point;

t = (dim + 1)*(dim + 2)/2;  %number of columns in matrix A

%Build A
temp = zeros(dim+1, size(M, 2));
for i = 1:size(M, 2)
    temp(1:dim, i) = M(i).Point;
    temp(dim+1, i) = M(i).Value;
end

A = buildA(temp);

%Define F
F = temp((dim + 1), 1:size(temp, 2))';

%Look for the minimizer of 1/2 *(Ay-F)' *(Ay-F) + tau/2 * (y' y)
%Rewrite the expression such that it fits into qudratic programing
% min(1/2 * y' * H * y + b' * y)
H = A' * A + tau * eye(t);
b = - A' * F;

if sum(isinf(H), 'all') + sum(isinf(b), 'all') == 0
    y = quadprog(H, b,[],[],[],[],[],[],[], optimset('Display', 'off'));
    
    %Transform y into a quadratic function with matrix Q and g
    Q = buildQ(y, dim);
    g = buildg(y, dim);
    
    %Find minimizer of x' * Q * x + g' * x
    %x = quadprog(Q, g,[],[],[],[],[],[],[], optimset('Display', 'off'));
    if sum(isinf(Q), 'all') + sum(isinf(g), 'all') + sum(isnan(Q), 'all') + sum(isnan(g), 'all') == 0
        x = trust(g,  Q, step * 2);
        if (isreal(x_Q.Point))
        x_Q.Point = x;
        grad = (Q + Q')*x + g;
        end
    end 
    
    
end
end