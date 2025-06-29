# Benchmarking N-body Jacobian


## N-body Jacobian with analytical method

```julia
@benchmark HighFidelityEphemerisModel.dfdx_Nbody_Interp([1.0, 0.0, 0.3, 0.5, 1.0, 0.0], 0.0, , 0.0)
```
```
BenchmarkTools.Trial: 10000 samples with 4 evaluations.
 Range (min … max):  7.552 μs …  15.104 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     7.636 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   7.684 μs ± 233.761 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▁▆██▇▇▆▄▃▂▁▁▁▂▁▁                                            ▂
  ██████████████████▇█▇▇▆▇▅▄▅▅▅▅▄▄▃▄▁▅▁▄▆██▇▇▆▆▄▁▄▄▅▅▆▁▁▁▅▁▄▅ █
  7.55 μs      Histogram: log(frequency) by time      8.88 μs <

 Memory estimate: 1.73 KiB, allocs estimate: 43.
```


## N-body Jacobian with ForwardDiff

```julia
@benchmark HighFidelityEphemerisModel.eom_jacobian_fd(HighFidelityEphemerisModel.eom_Nbody_Interp, [1.0, 0.0, 0.3, 0.5, 1.0, 0.0], 0.0, , 0.0)
```
```
BenchmarkTools.Trial: 10000 samples with 10 evaluations.
 Range (min … max):  1.913 μs … 480.854 μs  ┊ GC (min … max): 0.00% … 98.81%
 Time  (median):     2.158 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.464 μs ±   9.848 μs  ┊ GC (mean ± σ):  9.41% ±  2.42%

      ▁▅█▄▇▃▂                                                  
  ▁▁▃▃███████▇█▇▄▆▆▇▆▇▇▇▆██▅▆▆▅▄▄▃▃▃▃▃▂▃▃▃▂▃▃▂▂▃▂▂▂▂▂▁▂▁▂▁▁▁▁ ▃
  1.91 μs         Histogram: frequency by time        2.83 μs <

 Memory estimate: 6.69 KiB, allocs estimate: 50.
```

