function [V,a_policy,c_policy,info] = sub_vfi(params,numerics)
% SUB_VFI Solve the income fluctuation problem by value function iteration.
%
% INPUTS:
% params:   Structure containing model parameters, grids, and transition matrix
% numerics: Structure containing VFI and golden-search tolerances and limit
%
% OUTPUTS:
% V:        [n_a,n_z] value function
% a_policy: [n_a,n_z] next-period asset policy
% c_policy: [n_a,n_z] consumption policy
% info:     Structure containing convergence status, iterations, and error

a_grid = params.a_grid;
z_grid = params.z_grid;
pi_z = params.pi_z;
R = params.R;
beta = params.beta;
gamma = params.gamma;
n_a = params.n_a;
n_z = params.n_z;
a_min = params.a_min;
a_max = params.a_max;
b = params.b;
tol_vfi = numerics.tol_vfi;
max_iter = numerics.max_iter;
tol_golden = numerics.tol_golden;

V = zeros(n_a,n_z);
for iz = 1:n_z
    V(:,iz) = log(R*a_grid+z_grid(iz)+b)/(1-beta);
end
V_new = zeros(n_a,n_z);
a_policy = zeros(n_a,n_z);
c_floor = 1.0e-12;
error_v = Inf;
iteration = 0;

while error_v > tol_vfi && iteration < max_iter
    iteration = iteration + 1;

    % EV(ia,iz) = sum_izprime V(ia,izprime)*pi_z(iz,izprime).
    EV = V*pi_z.';

    for iz = 1:n_z
        EV_z = EV(:,iz);
        z_now = z_grid(iz);
        for ia = 1:n_a
            resources = R*a_grid(ia) + z_now;
            a_upper = min(a_max,resources-c_floor);

            % This closure is intentionally the benchmark whose JIT overhead
            % later implementations will seek to reduce.
            objective = @(aprime) rhs_bellman(aprime,resources,beta, ...
                gamma,a_grid,EV_z);
            [a_policy(ia,iz),V_new(ia,iz)] = golden(objective,a_min, ...
                a_upper,tol_golden);
        end
    end

    error_v = max(abs(V_new-V),[],'all');
    V(:) = V_new;
end

c_policy = R*a_grid + z_grid.' - a_policy;
info.converged = error_v <= tol_vfi;
info.iterations = iteration;
info.error = error_v;

end %end function
