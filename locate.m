function left = locate(x_grid,x_query)
% LOCATE Find the interval immediately to the left of a scalar query.
%
% INPUTS:
% x_grid:  [n,1] monotone increasing grid
% x_query: Scalar query point
%
% OUTPUTS:
% left:    Scalar interval index, with 0 or n indicating extrapolation

n = length(x_grid);
if x_query < x_grid(1)
    left = 0;
elseif x_query > x_grid(n)
    left = n;
else
    left = 1;
    right = n;
    while right-left > 1
        middle = floor((right+left)/2);
        if x_query >= x_grid(middle)
            left = middle;
        else
            right = middle;
        end
    end
end

end %end function

