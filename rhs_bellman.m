function value = rhs_bellman(aprime,resources,beta,gamma,a_grid,EV_z)
% RHS_BELLMAN Evaluate the Bellman objective for one scalar asset choice.
%
% INPUTS:
% aprime:    Scalar next-period asset choice
% resources: Scalar current resources, R*a + z
% beta:      Scalar discount factor
% gamma:     Scalar CRRA coefficient
% a_grid:    [n_a,1] asset grid
% EV_z:      [n_a,1] expected next-period value at grid points
%
% OUTPUTS:
% value:     Scalar Bellman objective value

consumption = resources-aprime;
if consumption <= 0
    value = -Inf;
    return
end

if gamma == 1
    utility = log(consumption);
else
    utility = consumption^(1-gamma)/(1-gamma);
end
value = utility + beta*interp1_scal(a_grid,EV_z,aprime);

end %end function

