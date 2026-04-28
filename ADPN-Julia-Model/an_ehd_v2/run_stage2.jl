# -----------------------------------------------------------------------------
# Stage 2 — activate D_mix(ε_org) sweep.
#   ε_org ∈ {0.02, 0.05, 0.09, 0.15, 0.25, 0.30}
#   δ   = 50 μm
# Spans single-phase (ε < ε_sat ≈ 0.0862) and two-phase (ε ≥ ε_sat).
# Exact Bloomquist experimental points + ε_org = 0.02 Stage 1 anchor.
# -----------------------------------------------------------------------------

include(joinpath(@__DIR__, "sweep_runner.jl"))

const STAGE2_DELTA_UM = 50.0
const STAGE2_EPS_ORG  = (0.02, 0.05, 0.09, 0.15, 0.25, 0.30)

println("\n#############################################################")
println("# STAGE 2: ε_org sweep at δ = $(STAGE2_DELTA_UM) μm")
println("# Values: ", STAGE2_EPS_ORG)
println("#############################################################\n")

for eps in STAGE2_EPS_ORG
    run_one_sweep(eps, STAGE2_DELTA_UM;
                  stage_tag = "stage2",
                  profile_V_targets = (-1.5, -1.8, -2.0, -2.5),
                  verbose = true)
    println("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")
end

println("\n\n#############################################################")
println("# Stage 2 complete.  Sweep of $(length(STAGE2_EPS_ORG)) ε_org values finished.")
println("#############################################################")
