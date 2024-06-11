# Small Path Tracer

## Features

- Monte Carlo Path Tracing
- Russian Roulette
- MultiSampling Anti-Aliasing (MSAA)
- Bounding Volume Hierarchy (BVH)
- Mircofacet Model
- Importance Sampling
- Multiple Importance Sampling (MIS)
- Multithreading

## Running
```
mkdir build
cd build
cmake -G "MinGW Makefiles" .. (first time on VS2022) / cmake ..
make 
./RayTracing
```

## Results

<center><img src="results/all.png"></center>


## Acknowledgement

The code is based on *[GAMES101](https://sites.cs.ucsb.edu/~lingqi/teaching/games101.html)*, *[GAMES202](https://sites.cs.ucsb.edu/~lingqi/teaching/games202.html)*, *[Ray Tracing in One Weekend Book Series](https://github.com/RayTracing/raytracing.github.io)*, *[Howl's Blog](https://howl144.github.io/2023/09/30/00014.%20Games101%20FinalProject/#shadowing-masking-function)* and *[YANGTHEKING](https://blog.csdn.net/ycrsw/article/details/124408789)*.
