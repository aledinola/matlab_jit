# MATLAB JIT benchmark: income fluctuation problem

This repository contains continuous-choice value function iteration (VFI) implementations of the income fluctuation problem in MATLAB. `sub_vfi` is the benchmark implementation. `sub_vfi_1` is the first improved implementation and inlines golden-section search inside the VFI state loop so MATLAB's JIT compiler sees the search control flow directly.

The benchmark intentionally creates an anonymous function of the form

```matlab
objective = @(aprime) rhs_bellman(aprime,resources,beta,gamma,a_grid,EV_z);
```

at each state and passes it to `golden`. This design makes the anonymous-function and function-handle overhead explicit for future comparisons.

`sub_vfi_1` preserves the same Bellman function, interpolation routine, search arithmetic, tolerances, and initialization. It replaces the call to `golden` and the anonymous objective with the inlined golden-section loop and direct calls such as

```matlab
f_left = rhs_bellman(x_left,resources,beta,gamma,a_grid,EV_z);
```

The golden-ratio constants, which are invariant across states and iterations, are also computed once per solver call in `sub_vfi_1`.

## Model and algorithm

The household solves

```text
V(a,z) = max_{a'} { u(R*a + z - a') + beta * E[V(a',z') | z] }
```

subject to positive consumption and asset-grid bounds, with log utility. Income follows a two-state Markov chain. Each VFI iteration precomputes conditional continuation values as

```matlab
EV = V*transpose(pi_z);
```

The maximization over continuous next-period assets uses golden-section search, and continuation values between asset grid points use scalar linear interpolation. There is no Howard acceleration.

The economic calibration follows QuantEcon's Julia lecture [Optimal Savings III: Occasionally Binding Constraints](https://julia.quantecon.org/dynamic_programming/ifp.html):

| Parameter | Value |
|---|---:|
| Gross return `R` | 1.01 |
| Discount factor `beta` | 0.96 |
| Utility | `log(c)` |
| Borrowing parameter `b` | 0 |
| Asset grid | 600 points on [0, 16] |
| Income states `z_grid` | [0.5, 1.0] |
| Transition matrix `pi_z` | [0.60 0.40; 0.05 0.95] |
| VFI tolerance | 1e-5 |
| Maximum VFI iterations | 10,000 |
| Golden-search tolerance | 1e-5 |

The Julia lecture's comparison uses 80 fixed VFI iterations rather than a convergence tolerance. This benchmark retains its economic calibration but uses a `1e-5` stopping tolerance and a generous 10,000-iteration safeguard so convergence can be tested explicitly.

## Run

From MATLAB in the repository directory:

```matlab
main
```

Run all utility and full-model tests with:

```matlab
run_tests
```

The full-model tests require convergence, a finite value function strictly increasing in assets, a next-period asset policy weakly increasing in current assets, strictly positive consumption, asset-bound feasibility, a row-stochastic transition matrix, and satisfaction of the budget identity.

## Measured running time

Measured on 2026-08-25 with MATLAB R2026a Update 4 (64-bit Windows) and the 600-point asset grid:

| Implementation | Warmed `timeit` runtime |
|---|---:|
| `sub_vfi` benchmark | **1.748796 seconds** |
| `sub_vfi_1` inline golden | **1.341655 seconds** |

Inlining golden-section search produced a measured **1.303x speedup**. Both implementations converged in 297 VFI iterations with final sup-norm error approximately `9.727e-6`, and their value and policy arrays agreed to the test tolerance of `1e-12`.

The measurements time only the complete solver calls; parameter construction, correctness tests, and console output are excluded. Timing is machine- and MATLAB-release-specific.

### JIT-aware benchmarking protocol

The reported `timeit` results measure steady-state execution after MATLAB has had an opportunity to JIT-compile both functions. They can be reproduced in one MATLAB session with:

```matlab
main  % Initialize inputs and validate both implementations

timing = benchmark_vfi(params,numerics);
```

`benchmark_vfi` first runs both implementations once as an untimed warm-up and then applies `timeit` to each. Both solvers initialize the value function on every invocation, so every timed evaluation performs a complete VFI solve. MATLAB reuses the compiled function code across evaluations. The creation and invocation of the anonymous Bellman objective inside `sub_vfi` are part of every benchmark solve and hence are included in its result. The outer zero-input wrappers required by `timeit` are included for both implementations.

The first-solve and warmed timings answer different questions:

- The first-solve time includes one-time loading and JIT-compilation effects.
- The warmed `timeit` result measures the recurring cost once compiled code is available and is the primary benchmark for studying JIT-sensitive implementations.

Do not run `clear functions` or `clear all` between warm-up and measurement, because either can discard compiled function state. The comparison uses identical inputs and measures both implementations in the same MATLAB session:

```matlab
% Warm up each complete solver once.
sub_vfi(params,numerics);
sub_vfi_1(params,numerics);

% Compare steady-state runtimes after JIT compilation.
t_benchmark = timeit(@() sub_vfi(params,numerics));
t_inline = timeit(@() sub_vfi_1(params,numerics));
```

Correctness should be checked separately before interpreting the timing comparison.

## Files

- `main.m`: calibration, grids, transition matrix, solve, timing, and quick tests.
- `sub_vfi.m`: benchmark VFI implementation.
- `sub_vfi_1.m`: improved VFI with inline golden-section search.
- `benchmark_vfi.m`: JIT warm-up and `timeit` comparison of both solvers.
- `rhs_bellman.m`: scalar Bellman objective.
- `golden.m`: scalar golden-section maximizer adapted from the MATLAB skill templates.
- `interp1_scal.m` and `locate.m`: scalar linear interpolation and interval search adapted from the MATLAB skill templates.
- `run_tests.m`: utility tests and full-model validation.
