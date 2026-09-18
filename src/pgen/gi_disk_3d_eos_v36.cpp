//========================================================================================
// Athena++ GI disk v3.6: v3.5 physics and mesh with optional inactive polar cones.
//========================================================================================
//! \file gi_disk_3d_eos_v36.cpp
//! \brief V3.5 disk with an opt-in theta_cut hydro mask for full-polar meshes.

#define GI_DISK_PGEN_NAME "gi_disk_3d_eos_v36"
#define GI_DISK_PRESSURE_CORRECTED_DEFAULT true
#define GI_DISK_ENABLE_MESHGEN true
#define GI_DISK_ENABLE_THETA_MASK true

#include "gi_disk_3d_eos_v2.cpp"
