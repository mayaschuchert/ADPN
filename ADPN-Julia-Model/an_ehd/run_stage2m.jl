# -----------------------------------------------------------------------------
# Stage 2m — m_i-corrected effective diffusivity (guide §4.2 upgrade).
#   Same ε_org sweep as Stage 2, but D_eff = (1−ε)·D_aq + ε·m_i·D_org
#   instead of the arithmetic mean. Only the two-phase regime changes;
#   single-phase still uses D_aq.
#
#   Partition coefficients (from params.jl m_partition):
#     m_AN = 11.59, m_ADPN = 7.72, m_PN = 11.59, ions = 0.
# -----------------------------------------------------------------------------

include(joinpath(@__DIR__, "sweep_runner.jl"))

const STAGE2M_DELTA_UM = 50.0
const STAGE2M_EPS_ORG  = (0.02, 0.05, 0.09, 0.15, 0.25, 0.30)

# Activate m-corrected formulation for this run.
set_D_formulation!(:m_corrected)

println("\n#############################################################")
println("# STAGE 2m: m_i-corrected D_eff, ε_org sweep at δ = $(STAGE2M_DELTA_UM) μm")
println("# D formulation: $(D_formulation())")
println("# Values: ", STAGE2M_EPS_ORG)
println("#############################################################\n")

for eps in STAGE2M_EPS_ORG
    run_one_sweep(eps, STAGE2M_DELTA_UM;
                  stage_tag = "stage2m",
                  profile_V_targets = (-1.5, -1.8, -2.0, -2.5),
                  verbose = true)
    println("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")
end

println("\n\n#############################################################")
println("# Stage 2m complete (m_i-corrected D_eff).")
println("# Compare with Stage 2 (arithmetic) using plot_stage2_comparison.py")
println("#############################################################")
