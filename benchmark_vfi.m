function timing = benchmark_vfi(params,numerics)
% BENCHMARK_VFI Compare JIT-warmed runtimes of the two VFI implementations.
%
% INPUTS:
% params:   Structure containing model parameters, grids, and transition matrix
% numerics: Structure containing VFI and golden-search tolerances and limit
%
% OUTPUTS:
% timing:   Structure containing both runtimes and benchmark/improved speedup

% Warm up both implementations before measuring steady-state execution.
sub_vfi(params,numerics);
sub_vfi_1(params,numerics);

timing.sub_vfi = timeit(@() sub_vfi(params,numerics));
timing.sub_vfi_1 = timeit(@() sub_vfi_1(params,numerics));
timing.speedup = timing.sub_vfi/timing.sub_vfi_1;

fprintf('sub_vfi:   %.9f seconds\n',timing.sub_vfi);
fprintf('sub_vfi_1: %.9f seconds\n',timing.sub_vfi_1);
fprintf('Speedup:   %.3fx\n',timing.speedup);

end %end function
