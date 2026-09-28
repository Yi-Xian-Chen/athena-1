//========================================================================================
// Athena++ GI disk v4.2: thermodynamically consistent floor states and an optional
// radius-dependent disk-monopole correction to the initial rotation profile.
//========================================================================================
//! \file gi_disk_3d_eos_v42.cpp
//! \brief V4 GI disk with Helmholtz-floor atmosphere and disk-monopole rotation.
// Recompile this wrapper whenever the shared v2 implementation changes.

#define GI_DISK_PGEN_NAME "gi_disk_3d_eos_v42"
#define GI_DISK_PRESSURE_CORRECTED_DEFAULT true
#define GI_DISK_ENABLE_MESHGEN true
#define GI_DISK_ENABLE_THETA_MASK true
#define GI_DISK_ATMOSPHERE_FROM_HELM_TFLOOR true

#include "gi_disk_3d_eos_v2.cpp"
