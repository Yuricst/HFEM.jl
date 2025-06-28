# Benchmarking N-body Jacobian


## N-body Jacobian with analytical method

```julia
@benchmark HFEM.dfdx_Nbody_Interp([1.0, 0.0, 0.3, 0.5, 1.0, 0.0], 0.0, , 0.0)
```
```
BenchmarkTools.Trial: 10000 samples with 4 evaluations.
 Range (min … max):  7.708 μs …  14.906 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     7.885 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   7.929 μs ± 256.359 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

   ▄█▂  ▁▄▅▆▄                                                  
  ▃███▇▇██████▆▄▃▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▁▁▁▂▁▁▁▁▁▁▁▁▁ ▂
  7.71 μs         Histogram: frequency by time           9 μs <

 Memory estimate: 1.73 KiB, allocs estimate: 43.
```


## N-body Jacobian with ForwardDiff

```julia
@benchmark HFEM.dfdx_Nbody_Interp_fd([1.0, 0.0, 0.3, 0.5, 1.0, 0.0], 0.0, , 0.0)
```
```
BenchmarkTools.Trial: 10000 samples with 9 evaluations.
 Range (min … max):  1.968 μs … 626.519 μs  ┊ GC (min … max): 0.00% … 99.10%
 Time  (median):     2.232 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.493 μs ±  11.566 μs  ┊ GC (mean ± σ):  9.22% ±  1.98%

      ▃▄█▄▃                                                    
  ▁▁▃▅█████▆▆▄▅▅▆▅▆▅▇█▇█▆▇▅▆▄▆▄▄▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▁▁▁▁▁▁▁▁ ▃
  1.97 μs         Histogram: frequency by time        2.92 μs <

 Memory estimate: 6.69 KiB, allocs estimate: 50.
```

