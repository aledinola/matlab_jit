%% Solve and validate the benchmark income fluctuation problem.
clear;
clc;

% Economic parameters: QuantEcon's Julia income fluctuation lecture.
params.r = 0.01;
params.R = 1+params.r;
params.beta = 0.96;
params.gamma = 1.0; % Log utility
params.b = 0.0;

% State grids and Markov transition matrix.
params.a_min = -params.b;
params.a_max = 16.0;
params.n_a = 50;
params.z_grid = [0.5;1.0];
params.pi_z = [0.60,0.40;0.05,0.95];
params.n_z = length(params.z_grid);
params.a_grid = linspace(params.a_min,params.a_max,params.n_a).';

% The Julia lecture uses 80 fixed iterations. A tolerance and generous cap
% are used here so the requested convergence test is meaningful.
numerics.tol_vfi = 1.0e-5;
numerics.max_iter = 10000;
numerics.tol_golden = 1.0e-5;

fprintf('Solving the benchmark income fluctuation problem...\n');
tic;
[V,a_policy,c_policy,info] = sub_vfi(params,numerics);
vfi_seconds = toc;

% Quick economic and numerical tests.
assert(info.converged,'VFI did not converge within the iteration limit.');
assert(all(isfinite(V),'all'),'The value function contains nonfinite values.');
assert(all(diff(V,1,1) > 0,'all'), ...
    'The value function is not strictly increasing in assets.');
assert(all(diff(a_policy,1,1) >= -10*numerics.tol_golden,'all'), ...
    'The asset policy is not weakly increasing in current assets.');
assert(all(c_policy > 0,'all'),'The consumption policy is not strictly positive.');
assert(all(a_policy >= params.a_min & a_policy <= params.a_max,'all'), ...
    'The asset policy violates the asset-grid bounds.');
assert(max(abs(sum(params.pi_z,2)-1)) < 1.0e-12, ...
    'Rows of the transition matrix do not sum to one.');
assert(all(params.pi_z >= 0,'all'), ...
    'The transition matrix contains negative probabilities.');

resources = params.R*params.a_grid + params.z_grid.';
budget_error = max(abs(c_policy-(resources-a_policy)),[],'all');
assert(budget_error < 1.0e-12,'The budget constraint is not satisfied.');

fprintf('Converged in %d iterations (sup-norm error %.3e).\n', ...
    info.iterations,info.error);
fprintf('VFI running time: %.6f seconds.\n',vfi_seconds);
fprintf('All tests passed.\n');
