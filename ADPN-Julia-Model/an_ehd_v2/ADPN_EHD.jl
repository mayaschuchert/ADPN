# Master module: load all submodules so `using .ADPN_EHD` (or include("ADPN_EHD.jl"))
# exposes the full API as a single namespace.
module ADPN_EHD

include("params.jl")
include("mesh.jl")
include("diffusivity.jl")
include("chemistry.jl")
include("kinetics.jl")
include("transport.jl")
include("assembly.jl")
include("solver.jl")

# v6 additions (§17–§20)
include("hydrodynamics.jl")
include("cell_voltage.jl")
include("fixed_j_solver.jl")
include("fit_kinetics.jl")

using .Params
using .Mesh
using .Diffusivity
using .Chemistry
using .Kinetics
using .Transport
using .Assembly
using .Solver
using .Hydrodynamics
using .CellVoltage
using .FixedJ
using .FitKinetics

# re-export a useful surface
export Params, Mesh, Diffusivity, Chemistry, Kinetics, Transport, Assembly, Solver,
       Hydrodynamics, CellVoltage, FixedJ, FitKinetics,
       make_mesh, D_mix, set_D_formulation!, D_formulation,
       solve_phosphate_equilibrium, buffer_sources!, make_initial_guess,
       bulk_concentration, c_AN_bulk,
       tafel_currents, sg_flux,
       set_kinetic_override!, clear_kinetic_override!, with_kinetic_override,
       conc_dof, phi_dof, full_residual!, decode_state,
       faradaic_currents_from_state,
       newton_solve!, newton_continuation,
       build_banded_pattern, banded_fd_jacobian!,
       # v6 §17 — cell-voltage decomposition
       kappa_dilute, kappa_eff, R_series,
       V_cell_predicted, V_cathode_target,
       V_CE_DEFAULT, R_CONTACT_DEFAULT,
       # v6 §18 — hydrodynamics
       ml_min_to_m3_s, d_hydraulic, v_super,
       reynolds, schmidt, sherwood_leveque,
       delta_leveque, weber_numbers,
       # v6 §20 — fixed-j solver and fit driver
       solve_at_j, FixedJResult,
       BloomquistRow, FitContext, build_context,
       select_core, select_extended, select_holdout,
       residuals!, loss, theta_to_physical, physical_to_theta,
       lm_fit, LMResult, N_THETA, THETA_LB, THETA_UB, THETA0

end # module
