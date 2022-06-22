%Author: Dominic (Zhongda) Huang
%Date: 2021.06.02
%Input: vector of coefficients y, dimension of input
%Output: vector g

function g = buildg(y, dim)
    g = y(end-dim:end-1);
end