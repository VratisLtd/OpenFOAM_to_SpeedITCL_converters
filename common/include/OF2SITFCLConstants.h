#ifndef OF_2_SITFCL_CONSTANTS_H
#define OF_2_SITFCL_CONSTANTS_H

namespace speeditcl
{
namespace io
{

////////////////////////////////////////////////////////////////////////////////

constexpr auto DATA_DIR_NAME      = "./sitfcl";
constexpr auto FILE_NAME_MESH     = "mesh";
constexpr auto FILE_NAME_SETTINGS = "settings";

////////////////////////////////////////////////////////////////////////////////

constexpr auto FILE_HEADER        = "SITFCL_FILE_1.0";

////////////////////////////////////////////////////////////////////////////////

constexpr auto VECTOR_HEADER      = "SITFCL_VECTOR";

constexpr auto DATA_DOUBLE        = "SITFCL_DOUBLE";
constexpr auto DATA_INT           = "SITFCL_INT";

////////////////////////////////////////////////////////////////////////////////

constexpr auto MODE_TEXT          = "SITFCL_TEXT_MODE";
constexpr auto MODE_BINARY        = "SITFCL_BINARY_MODE";

////////////////////////////////////////////////////////////////////////////////

constexpr auto FIELD_HEADER              = "SITFCL_FIELD";

constexpr auto CELL_VALS_HEADER          = "SITFCL_CELL_VALUES_HEADER";
constexpr auto BFACE_VALS_HEADER         = "SITFCL_BOUNDARY_FACE_VALUES_HEADER";
constexpr auto BFACE_TYPES_HEADER        = "SITFCL_BOUNDARY_FACE_TYPES_HEADER";
constexpr auto FIELD_OF_HEADER           = "SITFCL_OF_FIELD_HEADER";
constexpr auto FIELD_OF_FOOTER           = "SITFCL_OF_FIELD_FOOTER";
constexpr auto FIELD_OF_BOUNDARY_PATCHES = "SITFCL_OF_BOUNDARY_PATCHES";

////////////////////////////////////////////////////////////////////////////////

constexpr auto MESH_HEADER                 = "SITFCL_MESH";

constexpr auto IFACE_AREA_VECTORS_HEADER   = "SITFCL_IFACE_AREA_VECTORS";
constexpr auto IFACE_WEIGHTS_HEADER        = "SITFCL_IFACE_WEIGHTS";
constexpr auto IFACE_DELTA_COEFFS_HEADER   = "SITFCL_IFACE_DELTA_COEFFS";
constexpr auto IFACE_CORRECTION_VEC_HEADER = "SITFCL_IFACE_CORRECTION_VECTORS";
constexpr auto IFACE_OWNERS_HEADER         = "SITFCL_IFACE_OWNERS";
constexpr auto IFACE_NEIGHBOURS_HEADER     = "SITFCL_IFACE_NEIGHBOURS";

constexpr auto BFACE_AREA_VECTORS_HEADER   = "SITFCL_BFACE_AREA_VECTORS";
constexpr auto BFACE_WEIGHTS_HEADER        = "SITFCL_BFACE_WEIGHTS";
constexpr auto BFACE_DELTA_COEFFS_HEADER   = "SITFCL_BFACE_DELTA_COEFFS";
constexpr auto BFACE_CORRECTION_VEC_HEADER = "SITFCL_BFACE_CORRECTION_VECTORS";
constexpr auto BFACE_OWNERS_HEADER         = "SITFCL_BFACE_OWNERS";
constexpr auto BFACE_REORDER_TABLE_HEADER  = "SITFCL_BFACE_REORDER_TABLE";

constexpr auto CELL_VOLUMES_HEADER         = "SITFCL_CELL_VOLUMES";

////////////////////////////////////////////////////////////////////////////////

constexpr auto SETTINGS_HEADER             = "SITFCL_SETTINGS";

constexpr auto SETTING_SOLVER_NAME         = "solver_name";
constexpr auto SOLVER_ICO                 =  "ICO";
constexpr auto SOLVER_SIMPLE               = "SIMPLE";
constexpr auto SOLVER_PISO                 = "PISO";

constexpr auto SETTING_CASE_DIR            = "case_dir";
constexpr auto SETTING_N_TIME_STEPS        = "n_time_steps";
constexpr auto SETTING_DELTA_T             = "delta_t";
constexpr auto SETTING_ADJUSTABLE_RUN_TIME = "adjustable_run_time";
constexpr auto SETTING_WRITE_CONTROL       = "writeControl";
constexpr auto SETTING_WRITEINTERVAL       = "writeInterval";
constexpr auto SETTING_START_TIME          = "start_time";
constexpr auto SETTING_END_TIME            = "end_time";
constexpr auto SETTING_MAX_CO              = "max_Co";

constexpr auto SETTING_P_FIXED_CELL_IDX    = "p_fixed_cell_idx";
constexpr auto SETTING_P_FIXED_CELL_VALUE  = "p_fixed_cell_value";
constexpr auto SETTING_NU                  = "nu";

constexpr auto SETTING_TURBULENCE_MODE     = "turbulence_mode";
constexpr auto TURBULENCE_K_OMEGA_SST      = "kOmegaSST";
constexpr auto TURBULENCE_LAMINAR          = "laminar";

constexpr auto TURBULENCE_KOMEGASST_ALPHA_K1 = "alphaK1";
constexpr auto TURBULENCE_KOMEGASST_ALPHA_K2 = "alphaK2";

constexpr auto TURBULENCE_KOMEGASST_ALPHA_OMEGA1 = "alphaOmega1";
constexpr auto TURBULENCE_KOMEGASST_ALPHA_OMEGA2 = "alphaOmega2";

constexpr auto TURBULENCE_KOMEGASST_GAMMA1 = "gamma1";
constexpr auto TURBULENCE_KOMEGASST_GAMMA2 = "gamma2";

constexpr auto TURBULENCE_KOMEGASST_BETA1 = "beta1";
constexpr auto TURBULENCE_KOMEGASST_BETA2 = "beta2";
constexpr auto TURBULENCE_KOMEGASST_BETASTAR = "betaStar";

constexpr auto TURBULENCE_KOMEGASST_A1 = "a1";
constexpr auto TURBULENCE_KOMEGASST_C1 = "c1";
constexpr auto TURBULENCE_KOMEGASST_CMU = "Cmu";

constexpr auto TURBULENCE_KOMEGASST_KAPPA = "kappa";
constexpr auto TURBULENCE_KOMEGASST_E = "E";

constexpr auto TURBULENCE_KOMEGASST_YPLUS_LAM = "YPlusLam";


constexpr auto SETTING_N_NONORTH_CORR_ITER = "n_nonorthogonal_correction_iterations";
constexpr auto SETTING_N_PISO_CORR_ITER    = "n_PISO_correction_iterations";

constexpr auto SETTING_P_REQUIRED_RESIDUAL = "p_required_residual";
constexpr auto SETTING_U_REQUIRED_RESIDUAL = "U_required_residual";


constexpr auto SETTING_EQN_MAX_ITER       = "max_iterations";
constexpr auto SETTING_EQN_REQ_TOLERANCE  = "required_tolerance";
constexpr auto SETTING_EQN_REL_TOLERANCE  = "relative_tolerance";
constexpr auto SETTING_EQN_PRECOND        = "preconditioner";

//constexpr auto SETTING_PEQN_MAX_ITER       = "p_solver_max_iterations";
//constexpr auto SETTING_PEQN_REQ_TOLERANCE  = "p_solver_required_tolerance";
//constexpr auto SETTING_PEQN_REL_TOLERANCE  = "p_solver_relative_tolerance";
//constexpr auto SETTING_PEQN_PRECOND        = "p_solver_preconditioner";
constexpr auto SETTING_AMG_PRECOND_REBUILD_COUNTER = "amg_rebuild_counter";

constexpr auto SETTING_U_RELAXATION_FACTOR = "U_relaxation_factor";
constexpr auto SETTING_P_RELAXATION_FACTOR = "p_relaxation_factor";



////////////////////////////////////////////////////////////////////////////////

constexpr auto TIME_VARYING_BOUNDARY_HEADER             = "SITFCL_TIME_VARYING_BOUNDARY";
constexpr auto TIME_VARYING_BOUNDARY_PATCH_HEADER       = "SITFCL_TIME_VARYING_BOUNDARY_PATCH";
constexpr auto TIME_VARYING_PATCH_FACE_POSITIONS_HEADER = "SITFCL_TIME_VARYING_PATCH_FACES";
constexpr auto TIME_VARYING_PATCH_VALUES_HEADER         = "SITFCL_TIME_VARYING_PATCH_VALUES";
constexpr auto TIME_VARYING_TIMELINES_TIME_HEADER      = "SITFCL_TIME_VARYING_TIME_TIMELINES";

////////////////////////////////////////////////////////////////////////////////

constexpr auto GROOVY_BOUNDARY_HEADER             = "SITFCL_GROOVY_BOUNDARY";
constexpr auto GROOVY_BOUNDARY_PATCH_HEADER       = "SITFCL_GROOVY_BOUNDARY_PATCH";
constexpr auto GROOVY_PATCH_FACE_POSITIONS_HEADER = "SITFCL_GROOVY_PATCH_FACES";
constexpr auto GROOVY_PATCH_VALUES_HEADER         = "SITFCL_GROOVY_PATCH_VALUES";

////////////////////////////////////////////////////////////////////////////////

constexpr auto GROOVY_TIMELINES_HEADER            = "SITFCL_GROOVY_TIMELINES";
constexpr auto GROOVY_TIMELINES_PATCH_HEADER      = "SITFCL_GROOVY_PATCH_TIMELINES";
constexpr auto GROOVY_TIMELINES_TIME_HEADER       = "SITFCL_GROOVY_TIME_TIMELINES";
constexpr auto GROOVY_TIMELINES_VALUES_HEADER     = "SITFCL_GROOVY_VALUES_TIMELINES";

////////////////////////////////////////////////////////////////////////////////

constexpr auto TOTALPRESSURE_BOUNDARY_HEADER             = "SITFCL_TOTALPRESSURE_BOUNDARY";
constexpr auto TOTALPRESSURE_BOUNDARY_PATCH_HEADER       = "SITFCL_TOTALPRESSURE_BOUNDARY_PATCH";
constexpr auto TOTALPRESSURE_PATCH_FACE_POSITIONS_HEADER = "SITFCL_TOTALPRESSURE_PATCH_FACES";
constexpr auto TOTALPRESSURE_PATCH_VALUES_HEADER         = "SITFCL_TOTALPRESSURE_PATCH_VALUES";

////////////////////////////////////////////////////////////////////////////////

} // namespace io
} // namespace speeditcl

#endif // OF_2_SITFCL_CONSTANTS_H
