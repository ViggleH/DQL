%Author: Dominic (Zhongda) Huang
%Date: 2021.06.02
%Check if a given number n is even.

function ans = isEven(n)
    ans = rem(n, 2) == 0;
end