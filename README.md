# MATLAB JIT benchmark: income fluctuation problem

This repository contains a standard continuous-choice value function iteration (VFI) implementation of the income fluctuation problem in MATLAB. `sub_vfi` is the benchmark implementation of the VFI, and we plan to improve it in the future by testing formulations that allow MATLAB's JIT compiler to eliminate more function-call overhead.

The benchmark intentionally creates an anonymous function of the form

```matlab
objective = @(aprime) rhs_bellman(aprime,resources,beta,gamma,a_grid,EV_z);
```

at each state and passes it to `golden`. This design makes the anonymous-function and function-handle overhead explicit for future comparisons.

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
| Asset grid | 50 points on [0, 16] |
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

Measured on 2026-08-25 with MATLAB R2026a Update 4 (64-bit Windows):

- Warmed `timeit(@() sub_vfi(params,numerics))`: **0.113074 seconds**.
- The preceding validation solve, including the first call in that MATLAB session: **0.148306 seconds**.
- Convergence: **297 VFI iterations**, with final sup-norm error **9.721e-6**.

The benchmark times only `sub_vfi`; parameter construction, tests, and console output are excluded. Timing is machine- and MATLAB-release-specific.

### JIT-aware benchmarking protocol

The reported `timeit` result measures steady-state execution after MATLAB has had an opportunity to JIT-compile the functions. It was obtained in one MATLAB session with:

```matlab
main  % Complete first solve; initializes inputs and warms the JIT compiler

benchmark_seconds = timeit(@() sub_vfi(params,numerics));
```

`sub_vfi` initializes the value function on every invocation, so every `timeit` evaluation performs a complete VFI solve. MATLAB reuses the compiled function code across evaluations. The creation and invocation of the anonymous Bellman objective inside `sub_vfi` are part of every timed solve and hence are included in the result. The outer zero-input wrapper required by `timeit` is also included.

The first-solve and warmed timings answer different questions:

- The first-solve time includes one-time loading and JIT-compilation effects.
- The warmed `timeit` result measures the recurring cost once compiled code is available and is the primary benchmark for studying JIT-sensitive implementations.

Do not run `clear functions` or `clear all` between warm-up and measurement, because either can discard compiled function state. For a fair comparison with a future implementation, warm up both versions and measure them in the same MATLAB session using identical inputs:

```matlab
% Warm up each complete solver once.
sub_vfi(params,numerics);
sub_vfi_improved(params,numerics);

% Compare steady-state runtimes after JIT compilation.
t_benchmark = timeit(@() sub_vfi(params,numerics));
t_improved = timeit(@() sub_vfi_improved(params,numerics));
```

Correctness should be checked separately before interpreting the timing comparison.

## Files

- `main.m`: calibration, grids, transition matrix, solve, timing, and quick tests.
- `sub_vfi.m`: benchmark VFI implementation.
- `rhs_bellman.m`: scalar Bellman objective.
- `golden.m`: scalar golden-section maximizer adapted from the MATLAB skill templates.
- `interp1_scal.m` and `locate.m`: scalar linear interpolation and interval search adapted from the MATLAB skill templates.
- `run_tests.m`: utility tests and full-model validation.
