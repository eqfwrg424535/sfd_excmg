SFD_EXCMG

Author: Dr. Jinxuan Wang

Affliation: China University of Geoscience

Release date: 10/12/24

Introduction

The SFD_EXCMG is an open-source program for large-scale MT forward modeling, which utilizes an extrapolation multigrid method to accelerate the solving of linear systems arising from staggered-grid 
finite difference (SFD) discretization of the curl-curl equation. The program is developed for complex geo-electrical settings, thus it supports arbitrary anisotropic conductivity and undulating 
topography. One can control the total number of levels and the tolerance on the finest mesh to see the differences. The formulations for calculating the background field are from Josef and Fernando(2002).
One can refer to the paper by Pan et al.(2024) for more details of the novel extrapolation multigrid method.

Installation

The SFD_EXCMG code is written with Fortran 90, and one can easily build the program on Windows or Linux devices.

For Windows users, they can build the program using powerful tools and softwares like Visual Studio, equipped with some certain Fortran compilers. We recommend the users to use Intel compilers like Intel Visual Fortran Compiler XE 19.0 or OneAPI Toolkit (2021.0 or newer), since the program is dependent on the direct solver mkl_pardiso from Intel libraries.

For Linux users, the program can be built using standard Makefile and make command. And again, Intel compilers should be installed in advance. The Makefile is already contained in the package. Before making, one should modify the following paths to their own, in case of compilation error.

mkllib=/opt/intel/oneapi/mkl/2022.0.2/lib/intel64

mklinc=/opt/intel/oneapi/mkl/2022.0.2/include

Please double check that you have compiled the code; that you are running the code on the same system as the one you used to compile it; that you are running the code correctly; that your file formats are correct. Please refer to the Documentation for details. 

If you have any questions, please feel free to contact us with the following e-mails:

jinxuanwang@cug.edu.cn/ jinx1219741@163.com 

References

Josef, P., Fernando A.M., S., 2002. Magnetotelluric impedances and parametric sensitivities for 1-D anisotropic layered media. Computers & Geosciences 28, 939–950.409.

Pan K., Wang J., Liu Z., Ou Z., Guo R., Ren Z. 2024. An efficient cascadic multigrid method combined with regularization technique for 3-D electromagnetic finite-element anisotropic modeling. Geophysics, 89(6), E241-E253.
