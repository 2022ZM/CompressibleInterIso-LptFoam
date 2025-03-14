# compressibleInterIsoLptFoam

## Features

- Volume-of-Fluid (VoF) solver for compressible flows coupled with Lagrangian library
- Conversion of small VoF elements to Lagrangian parcels to reduce computational cost
- Support for
  - Particle-particle interaction (collision, coalescence) and secondary breakup
  - Adaptive mesh refinement
  - Fully parallelized
- Based on interIsoFoam in [OpenFOAM v2112](https://www.openfoam.com).

## Compilation

 1. Clone the `CompressibleInterIso-LptFoam` gitHub repository:
```
git clone https://github.com/2022ZM/CompressibleInterIso-LptFoam.git
```
 2. Use the `Allwmake` script to compile the `compressibleInterIsoLptFoam` solver and its corresponding library:
```
./Allwmake
```

## Reference

The source code has been published in the following open-access research article:
```
Reference1
Meng Zhang, Xiang Lv, Fangmian Dong, et al. 
High-Fidelity atomization simulation of kerosene swirl injector using Multi-Resolution framework
Physics of Fluids，March 2025.
https://doi.org/10.1063/5.0258236

```
