function [x_best,f_best] = golden(objective,a_lower,a_upper,tolerance)
% GOLDEN Maximize a scalar objective over a bounded interval.
%
% INPUTS:
% objective: Function handle accepting one scalar and returning one scalar
% a_lower:   Scalar lower endpoint
% a_upper:   Scalar upper endpoint
% tolerance: Scalar absolute stopping tolerance for the search interval
%
% OUTPUTS:
% x_best:    Scalar approximate maximizer
% f_best:    Scalar objective value at x_best

alpha_left = (3-sqrt(5))/2;
alpha_right = (sqrt(5)-1)/2;
distance = a_upper-a_lower;
x_left = a_lower+alpha_left*distance;
x_right = a_lower+alpha_right*distance;
f_left = objective(x_left);
f_right = objective(x_right);

distance = alpha_left*alpha_right*distance;
while distance > tolerance
    distance = distance*alpha_right;
    if f_right < f_left
        x_right = x_left;
        x_left = x_left-distance;
        f_right = f_left;
        f_left = objective(x_left);
    else
        x_left = x_right;
        x_right = x_right+distance;
        f_left = f_right;
        f_right = objective(x_right);
    end
end

if f_right > f_left
    x_best = x_right;
    f_best = f_right;
else
    x_best = x_left;
    f_best = f_left;
end

end %end function

