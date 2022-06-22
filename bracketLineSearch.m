%Author: Dominic (Zhongda) Huang
%Date: 2021.06.29
%Input: Objective function f, three points that forms a bracket, number of
%iternation n.
%Output: best solution by safeguarded bracket linear search, x_Line and f_Line; the evaluated
%point-value pair list XF_Line

function [x_Line, f_Line, XF_Line] = bracketLineSearch(F, a_l, F_l, a_0, F_0, a_r, F_r, n)

%safe guarded parameter s.
s = 1/16 * min([1/(4*(a_r - a_l)), -a_l/(a_r-a_l)^2, a_r/(a_r-a_l)^2]);

XF_Line = zeros(2,n);
x_Line = [];
f_Line = [];

for i = 1:n
    %Checking bracket length
    if(abs(a_r - a_l) < 10^(-3))
        break;
    end
    
    
    p = polyfit([a_l, a_0, a_r], [F_l, F_0, F_r], 2);
    
    
    %Checking polyfit result
    if(sum(isnan(p)) + sum(isinf(p)) ~= 0)
        break;
    end
    
    if(p(1) == 0)
        R = a_0;
    else
        R = a_0 - (p(2) + 2*p(1)*a_0)/2*p(1);
    end
    sig = s*(a_r - a_l)^2;
    a = proj(R, [a_l + sig, a_r - sig]);
    if norm(a_0 - a) < sig
        if a_0 <= 1/2 * (a_l + a_r)
            a = a_0 + sig;
        else
            a = a_0 - sig;
        end
    end
    F_a = F(a);
    XF_Line(1,i) = a;
    XF_Line(2,i) = F_a;
    
    if(a > a_0)
        if(F_a >= F_0)
            x_Line = a_0;
            f_Line = F_0;
            
            a_r = a;
            F_r = F_a;
            
        else
            x_Line = a;
            f_Line = F_a;
            
            a_l = a_0; F_l = F_0;
            a_0 = a; F_0 = F_a;
            
        end
    else
        if(F_a >= F_0)
            x_Line = a_0;
            f_Line = F_0;
            
            a_l = a;
            F_l = F_a;
        else
            x_Line = a;
            f_Line = F_a;
            
            a_r = a_0; F_r = F_0;
            a_0 = a; F_0 = F_a;
        end
        
    end
    
end

%Trim unused entries
XF_Line = XF_Line(:,1:i);

end