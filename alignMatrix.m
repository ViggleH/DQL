%Author: Dominic (Zhongda) Huang
%Date: 2021.06.02
%Input: origianl matrix D1, to-be-aligned direction d2 
%Output: a rotation matrix that rotates D1 such that the first vector align
%with direction d2.

function M = alignMatrix(D1, d2)
    
    d1 = D1(:,1);
    if ((norm(d1) * norm(d2)) ~= 0)
        theta = acos(dot(d1, d2) / (norm(d1) * norm(d2)));
    else
        theta = 0;
    end
    
    if(theta == 0)
        M = eye(size(d1, 1));
    else
        if(theta == pi)
            M = -eye(size(d1, 1));
        else
        M = rotMatrix(theta, d1, d2);
        end
    end
end