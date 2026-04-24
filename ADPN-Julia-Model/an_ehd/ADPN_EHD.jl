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

using .Params
using .Mesh
using .Diffusivity
using .Chemistry
using .Kinetics
using .Transport
using .Assembly
using .Solver

# re-export a useful surface
export Params, Mesh, Diffusivity, Chemistry, Kinetics, Transport, Assembly, Solver,
       make_mesh, D_mix, set_D_formulation!, D_formulation,
       solve_phosphate_equilibrium, buffer_sources!, make_initial_guess,
       bulk_concentration, c_AN_bulk,
       tafel_currents, sg_flux,
       conc_dof, phi_dof, full_residual!, decode_state,
       faradaic_currents_from_state,
       newton_solve!, newton_continuation,
       build_banded_pattern, banded_fd_jacobian!

end # module
