%Author: Dominic (Zhongda) Huang
%Date: 2021.06.02
%Input: list of evaluated point-value pairs M
%Output: best solution by quadratic step, x_Q and f_Q; the evaluated
%point-value pair list XF_Q

function A = buildA(M)

    %Initialization
    dim = size(M, 1) - 1;  %dimension of the problem
    t = (dim + 1)*(dim + 2)/2;  %number of rows in matrix A
    A = ones(size(M,2), t);
    
    for i = 1:size(M,2)
        %write square terms
        for j = 1:dim
            A(i,j) = M(j ,i)^2;
        end
        
        %write cross terms
        if(dim >= 2)
            combos = nchoosek(M(1:dim,i),2);
            for j = 1:size(combos, 1)
                A(i, j + dim) = 2*combos(j,1)*combos(j,2);
            end
        else
            combos = M(1, i);
            A(i, j + dim) = 2 * combos(1);
        end
        
        
        %write constant terms
        for j = 1:dim
            A(i,j + dim + size(combos, 1)) = M(j, i);
        end
    end
    
end