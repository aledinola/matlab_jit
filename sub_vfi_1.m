function [V,a_policy,c_policy,info] = sub_vfi_1(params,numerics)
% SUB_VFI_1 Solve the income fluctuation problem with inline golden search.
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

% Golden-section constants are invariant across all state optimizations.
alpha_left = (3-sqrt(5))/2;
alpha_right = (sqrt(5)-1)/2;
alpha_product = alpha_left*alpha_right;

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

            % Inline golden-section maximization. The Bellman RHS remains a
            % separate function, but no function handle or closure is used.
            distance = a_upper-a_min;
            x_left = a_min+alpha_left*distance;
            x_right = a_min+alpha_right*distance;
            f_left = rhs_bellman(x_left,resources,beta,gamma,a_grid,EV_z);
            f_right = rhs_bellman(x_right,resources,beta,gamma,a_grid,EV_z);

            distance = alpha_product*distance;
            while distance > tol_golden
                distance = distance*alpha_right;
                if f_right < f_left
                    x_right = x_left;
                    x_left = x_left-distance;
                    f_right = f_left;
                    f_left = rhs_bellman(x_left,resources,beta,gamma, ...
                        a_grid,EV_z);
                else
                    x_left = x_right;
                    x_right = x_right+distance;
                    f_left = f_right;
                    f_right = rhs_bellman(x_right,resources,beta,gamma, ...
                        a_grid,EV_z);
                end
            end

            if f_right > f_left
                a_policy(ia,iz) = x_right;
                V_new(ia,iz) = f_right;
            else
                a_policy(ia,iz) = x_left;
                V_new(ia,iz) = f_left;
            end
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
