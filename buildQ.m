%Author: Dominic (Zhongda) Huang
%Date: 2021.06.02
%Input: vector of coefficients y, dimension of input
%Output: matrix Q


function Q = buildQ(y, dim)
    %diag terms
    Q = diag(y(1:dim));
    
    k = dim+1;
    for i = 1:dim-1
        for j = i+1:dim
            Q(i,j) = y(k);
            k = k + 1;
        end
    end
    Q = 2*((Q + Q') - eye(size(Q,1)).*diag(Q));
end