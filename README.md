# PASC_inference

HPC C++ library for signal denoising and clustering time-series. 
Based on FEM-H1 regularisation proposed by Illia Horenko et al. and supported by "PASC project: Towards the HPC-inference of causality networks from multiscale economical data".

Generated Doxyden documentation can be found here: http://eth-cscs.github.io/PASC_inference/

## Installation
- download library
```
git clone https://github.com/eth-cscs/PASC_inference.git
cd PASC_inference
```
- download submodules
```
./util/update_submodules
```
- be sure that right modules are loaded and PETSc library is installed (and `PETSC_DIR`, `PETSC_ARCH` are properly set); if you are compiling the code on PIZ Daint, you can use scripts, which automatically set right modules and PETSc
```
source util/module_load_daint
```
- compile examples
```
cd examples
mkdir build
cd build
cmake -DFIND_PETSC=ON ..
```
(in case of CUDA implementation, use `-DUSE_CUDA=ON`)
- see the list of avaiable examples and choose which one to compile using `-DTEST_...=ON`
```
cmake -DTEST_SIGNAL1D=ON ..
```
- compile the library with examples
```
make
```
