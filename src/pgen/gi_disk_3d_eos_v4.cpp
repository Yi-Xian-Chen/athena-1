//========================================================================================
// Athena++ GI disk v4: v3.5 disk physics and polar mesh with spherical self-gravity.
// Configure with --grav=sph -mpi -fft and an Eigen include path.
//========================================================================================
//! \file gi_disk_3d_eos_v4.cpp
//! \brief Alpha-Helmholtz GI disk with nonuniform full-polar mesh and self-gravity.

#define GI_DISK_PGEN_NAME "gi_disk_3d_eos_v4"
#define GI_DISK_PRESSURE_CORRECTED_DEFAULT true
#define GI_DISK_ENABLE_MESHGEN true

#include "gi_disk_3d_eos_v2.cpp"
