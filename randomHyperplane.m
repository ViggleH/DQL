%Author: Dominic (Zhongda) Huang
%Date: 2021.06.02
%Input: dimension of the hyperplane
%Output: random hyperplane defined by the vector a and b


function [a, b] = randomHyperplane(dim)
    
    %Initialization
    a = zeros(dim, 1);
    b = zeros(dim ,1);
    r = randperm(dim);
    
    a(r(1)) = rand(1);
    b(r(2)) = rand(1);
end