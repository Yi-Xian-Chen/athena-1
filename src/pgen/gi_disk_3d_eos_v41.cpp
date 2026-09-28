//========================================================================================
// Athena++ GI disk v4.1: v4 self-gravitating disk with a configurable,
// thermodynamically consistent Helmholtz atmosphere temperature floor.
//========================================================================================
//! \file gi_disk_3d_eos_v41.cpp
//! \brief V4 GI disk with atmosphere and mask states tied to the Helmholtz T floor.
// Recompile this wrapper whenever the shared v2 implementation changes.

#define GI_DISK_PGEN_NAME "gi_disk_3d_eos_v41"
#define GI_DISK_PRESSURE_CORRECTED_DEFAULT true
#define GI_DISK_ENABLE_MESHGEN true
#define GI_DISK_ENABLE_THETA_MASK true
#define GI_DISK_ATMOSPHERE_FROM_HELM_TFLOOR true

#include "gi_disk_3d_eos_v2.cpp"
