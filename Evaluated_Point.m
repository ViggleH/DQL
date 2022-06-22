%Author: Dominic (Zhongda) Huang
%Date: 2021.08.09
%Object for evaluated points
%Properties:
%   - Point: To-be-evaluated point
%   - Value: calculated value
%   - Gradient: The centered simplex gradient of current point
%   - Type: 0 - only value is known
%           1 - both value and gradient is known
%           2 - gradient and value


classdef Evaluated_Point
    properties
        Point
        Value
        Gradient
        Hessian
        Type
    end
    methods
        %Evaluate the function value at current point
        function obj = eval(f, obj)
            obj.Value = f(obj.Point);
            obj.Type = 0;
        end
        %Approximate the gradient at current point
        function obj = grad(obj, obj_list)
            dim = size(obj.Point, 1); %dimension of the point
            num = size(obj_list, 2); %number of the points used for approximation
            T = zeros(dim, num);   %Initiate direction matrix
            d_s = zeros(num, 1); %Initiate delta_s matrix
            for i = 1:num
                T(:,i) = obj_list(i).Point - obj.Point;
                d_s(i) = obj_list(i).Value - obj.Value;
            end
            if sum(isnan(T), 'all') + sum(isinf(T), 'all') == 0
                grad = pinv(T')*d_s;
                if sum(isnan(grad), 'all') + sum(isinf(grad), 'all') == 0
                    obj.Gradient = grad;
                    obj.Type = 1;
                end
            end
        end
        %Approximate the Hessian at current point
        function obj = hes(obj, obj_list)
            dim = size(obj.Point, 1); %dimension of the point
            num = size(obj_list, 2); %number of the points used for approximation
            T = zeros(dim, num);   %Initiate direction matrix
            d_s = zeros(num, dim); %Initiate delta_s matrix
            for i = 1:num
                T(:,i) = obj_list(i).Point - obj.Point;
                d_s(i,:) = (obj_list(i).Gradient - obj.Gradient)';
            end
            if sum(isnan(T), 'all') + sum(isinf(T), 'all') == 0
            hes = pinv(T')*d_s;
            if sum(isnan(hes), 'all') + sum(isinf(hes), 'all') == 0
                obj.Hessian = hes;
                obj.Type = 2;
            end
            end
        end
    end
end