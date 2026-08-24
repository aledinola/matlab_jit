function y_query = interp1_scal(x_grid,y_grid,x_query)
% INTERP1_SCAL Linearly interpolate one scalar query on a monotone grid.
%
% INPUTS:
% x_grid:  [n,1] strictly increasing interpolation grid
% y_grid:  [n,1] values at interpolation grid points
% x_query: Scalar query point
%
% OUTPUTS:
% y_query: Scalar linearly interpolated value

n = size(x_grid,1);
left = locate(x_grid,x_query);
left = max(min(left,n-1),1);
slope = (y_grid(left+1)-y_grid(left)) / ...
    (x_grid(left+1)-x_grid(left));
y_query = y_grid(left)+(x_query-x_grid(left))*slope;

end %end function

