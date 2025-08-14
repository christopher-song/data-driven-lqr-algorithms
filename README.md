# Data-driven LQR algorithms with interval and sampled data
This repository contains the code to generate the examples from the paper "[All Data-Driven LQR Algorithms Require at Least as Much Interval Data as System Identification](https://ieeexplore.ieee.org/document/11071970)". The code implements the offline ADP and system identification algorithms for three example systems with interval and sampled data. These algorithms are detailed in the aforementioned paper, and a complete description of the offline learning ADP algorithm may be found in Chapter 2 of the book "[Robust Adaptive Dynamic Programming](https://yu-jiang.github.io/radpbook/)".
## Structure
The main script for interval data examples is `main_intervals.m`. This script generates data from a variable number of intervals and runs ADP and system identification for a particular system to compute the optimal $K$ and $P$ matrices using the `adp_intervals` and `sysid_intervals` methods respectively. 

The main script for sampled data examples is `main_samples.m`.  This script runs ADP and system identification for a particular system to compute the optimal $K$ and $P$ matrices, with a variable number of samples per interval, with a variable number of intervals, with trapezoid and chebyshev quadratures. It does this using `adp_samples` and `sysid_samples` and generates the data these methods require using `getdata_trapezoid` and `getdata_chebyshev`. 

The implementations of the algorithms themselves are found in `adp.m` and `sysid.m`. Example 2-dimensional, 3-dimensional, and 6-dimensional systems are included. 

## Requirements

This code was written with MATLAB R2022a. The sampled data examples use the [chebfun](https://www.chebfun.org/download/) package.

## Results

The plots for the 3-dimensional system are shown below.

Interval data:
<img width="1067" height="800" alt="3d" src="https://github.com/user-attachments/assets/89cf28eb-c54d-45c2-99a7-ffc531683c03" />
Sampled data:
<img width="1067" height="800" alt="3ds1" src="https://github.com/user-attachments/assets/b34841e0-5512-4c60-bb0f-dbb4c56fafd4" />
