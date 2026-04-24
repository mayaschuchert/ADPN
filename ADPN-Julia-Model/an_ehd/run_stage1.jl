# -----------------------------------------------------------------------------
# Stage 1 — single-phase reference at ε_org = 0.02, δ = 50 μm.
# -----------------------------------------------------------------------------
# Thin wrapper around sweep_runner.jl. See that file for the pipeline logic.

include(joinpath(@__DIR__, "sweep_runner.jl"))

const STAGE1_EPS      = 0.02       # minimum-physical single-phase loading
const STAGE1_DELTA_UM = 50.0

run_one_sweep(STAGE1_EPS, STAGE1_DELTA_UM;
              stage_tag = "stage1",
              profile_V_targets = (-1.5, -1.8, -2.0, -2.5),
              verbose = true)
