# -----------------------------------------------------------------------------
# run_stage3.jl — pre-compute warm-start cache for Stage 4.
#
# Loops over the unique (gap, Q_total, ε_org) keys derived from the Bloomquist
# dataset, computes δ via Lévêque, runs the v5 bootstrap pipeline (α_buf and
# α_kin ramps at V_warm = -1.0 V), and writes the converged DOF vector to disk
# in the v5-compatible cache filename format:
#
#     output/cache/s_eo<ε>_d<δ_μm>_V-1.000000.bin
#
# After Stage 3, every Stage 4 LM iteration cold-starts at the same point as
# the warmstarts and avoids the V-walk for the first iteration. Persisted cache
# survives Julia restarts, so re-running run_stage4.jl from a fresh process
# stays fast.
#
# Skips keys whose cache file already exists. Safe to re-run any time.
# -----------------------------------------------------------------------------
using Printf
using Dates
using DelimitedFiles

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.Solver: newton_solve!

const DATA_FILE = joinpath(@__DIR__, "..", "Experimental_data", "bloomquist_data.csv")
const CACHE_DIR = joinpath(@__DIR__, "output", "cache")
const STAGE3_DIR = joinpath(@__DIR__, "output", "stage3")
const STAGE3_DATA = joinpath(STAGE3_DIR, "data")
const STAGE3_LOGS = joinpath(STAGE3_DIR, "logs")
isdir(CACHE_DIR)   || mkpath(CACHE_DIR)
isdir(STAGE3_DATA) || mkpath(STAGE3_DATA)
isdir(STAGE3_LOGS) || mkpath(STAGE3_LOGS)

const N_MESH  = 100
const STRETCH = 10.0
const V_WARM  = -1.0   # bootstrap target V vs SHE — same as sweep_runner

# ---------- Inlined bootstrap (avoids include of sweep_runner.jl, which re-includes
# ADPN_EHD and breaks the module binding chain in Main scope). ----------
function _residual_closure(mesh, eps_org, V, alpha_buf, alpha_kin, c_eq)
    return (F, u) -> full_residual!(F, u, mesh, eps_org, V,
                                    alpha_buf, alpha_kin, c_eq)
end

function bootstrap!(u::Vector{Float64}, mesh, eps_org::Float64, c_eq;
                    verbose::Bool = false)
    # Step 1 — verify zero-residual initial state at α_buf = α_kin = 0
    res! = _residual_closure(mesh, eps_org, V_WARM, 0.0, 0.0, c_eq)
    F0 = zeros(length(u))
    res!(F0, u)
    @assert maximum(abs.(F0)) < 1e-8 "Initial residual not ~0"

    # Step 2 — buffer ramp (α_buf 0 → 1 in 10 uniform steps)
    for α in range(0.1, 1.0; length = 10)
        res! = _residual_closure(mesh, eps_org, V_WARM, α, 0.0, c_eq)
        r = newton_solve!(u, res!; max_iter = 25, verbose = false)
        @assert r.converged "buffer ramp failed at α_buf=$α (|F|=$(r.normF))"
    end

    # Step 3 — kinetics ramp (α_kin geometric ×2 from 1e-6 → 1.0, 21 steps)
    α_list = [min(1.0e-6 * 2.0^k, 1.0) for k in 0:20]
    α_list[end] = 1.0
    for α in α_list
        res! = _residual_closure(mesh, eps_org, V_WARM, 1.0, α, c_eq)
        r = newton_solve!(u, res!; max_iter = 25, verbose = false)
        @assert r.converged "kinetics ramp failed at α_kin=$α (|F|=$(r.normF))"
    end
    return u
end

function save_solution(u::Vector{Float64}, eps_org::Float64,
                       delta_um::Real, V::Float64)
    fn = @sprintf("s_eo%.3f_d%d_V%.6f.bin", eps_org, round(Int, delta_um), V)
    path = joinpath(CACHE_DIR, fn)
    open(path, "w") do f
        write(f, Int64(length(u)))
        write(f, u)
    end
    return path
end

# ---------- Filename helper (v5-compatible) ----------
function cache_filename(eps_org::Float64, delta_m::Float64, V::Float64)
    delta_um = round(Int, delta_m * 1e6)
    return @sprintf("s_eo%.3f_d%d_V%.6f.bin", eps_org, delta_um, V)
end

cache_path(eps_org, delta_m, V) =
    joinpath(@__DIR__, "output", "cache", cache_filename(eps_org, delta_m, V))

# ---------- Build the unique (gap, Q_total, ε_org) key list ----------
function load_keys(data_csv::String; eps_min::Float64 = 0.04)
    raw, hdr = readdlm(data_csv, ','; header = true)
    hdr = vec(hdr)
    col(name) = findfirst(==(string(name)), hdr)
    keys_seen = Set{Tuple{Float64,Float64,Float64}}()
    keys = Tuple{Float64,Float64,Float64}[]
    for i in 1:size(raw, 1)
        gap_mm  = Float64(raw[i, col("gap_mm")])
        Q_total = Float64(raw[i, col("Q_total_mL_min")])
        eps_org = Float64(raw[i, col("phi_AN")])
        eps_org < eps_min && continue
        key = (gap_mm, Q_total, eps_org)
        key in keys_seen && continue
        push!(keys_seen, key)
        push!(keys, key)
    end
    return keys
end

# ---------- Bootstrap one key ----------
function build_warmstart(gap_mm::Float64, Q_total_mL::Float64, eps_org::Float64,
                          c_eq; verbose::Bool = false)
    gap_m = gap_mm * 1.0e-3
    delta_m = delta_leveque(gap_m, ml_min_to_m3_s(Q_total_mL))
    cache_fp = cache_path(eps_org, delta_m, V_WARM)

    if isfile(cache_fp)
        return (skipped = true, delta_m = delta_m, path = cache_fp,
                converged = true, note = "already cached")
    end

    mesh = make_mesh(N_MESH, delta_m; stretch = STRETCH)
    u = make_initial_guess(N_MESH, c_eq, eps_org)

    t0 = time()
    bootstrap!(u, mesh, eps_org, c_eq; verbose = verbose)

    # bootstrap! converges at V_BOOT = V_WARM with α_buf = α_kin = 1.
    # Verify the final state is at the expected V.
    res! = (F, x) -> full_residual!(F, x, mesh, eps_org, V_WARM, 1.0, 1.0, c_eq)
    F0 = zeros(10 * N_MESH)
    res!(F0, u)
    normF = maximum(abs.(F0))

    save_solution(u, eps_org, delta_m * 1e6, V_WARM)
    dt = time() - t0
    return (skipped = false, delta_m = delta_m, path = cache_fp,
            converged = normF < 1e-4,
            normF = normF, wall_s = dt,
            note = normF < 1e-4 ? "" : "post-bootstrap |F|∞=$(normF)")
end

# ---------- Main ----------
function main()
    println("=" ^ 72)
    println(" Stage 3 — warm-start cache builder for Stage 4 (v6 §12)")
    println(" $(now())")
    println("=" ^ 72)

    keys = load_keys(DATA_FILE)
    @printf("Loaded %d unique (gap, Q_total, ε_org) keys (ε_org ≥ 0.04)\n",
            length(keys))

    c_eq = solve_phosphate_equilibrium()
    @printf("Bulk equilibrium: pH = %.3f, c_OH = %.2f, c_PO4 = %.2f mol/m³\n\n",
            c_eq.pH, c_eq.OH, c_eq.PO4)

    index_path = joinpath(STAGE3_DATA, "stage3_warmstart_index.csv")
    open(index_path, "w") do io
        println(io, "gap_mm,Q_total_mL_min,phi_AN,delta_um,V_warm_V," *
                    "cache_file,converged,wall_s,note")
        n_built = 0
        n_skip  = 0
        n_fail  = 0
        for (i, (gap_mm, Q_total, eps_org)) in enumerate(keys)
            @printf("[%2d/%d] gap=%.2f mm  Q=%.0f mL/min  ε=%.4f ... ",
                    i, length(keys), gap_mm, Q_total, eps_org)
            try
                r = build_warmstart(gap_mm, Q_total, eps_org, c_eq;
                                    verbose = false)
                if r.skipped
                    @printf("SKIP (already cached, δ=%.1f μm)\n", r.delta_m * 1e6)
                    n_skip += 1
                    @printf(io, "%.2f,%.0f,%.4f,%.1f,%.6f,%s,%s,%.2f,%s\n",
                            gap_mm, Q_total, eps_org, r.delta_m * 1e6,
                            V_WARM, basename(r.path), r.converged, 0.0, r.note)
                elseif r.converged
                    @printf("OK   δ=%.1f μm  |F|∞=%.2e  %.1f s\n",
                            r.delta_m * 1e6, r.normF, r.wall_s)
                    n_built += 1
                    @printf(io, "%.2f,%.0f,%.4f,%.1f,%.6f,%s,%s,%.2f,%s\n",
                            gap_mm, Q_total, eps_org, r.delta_m * 1e6,
                            V_WARM, basename(r.path), true, r.wall_s, r.note)
                else
                    @printf("FAIL |F|∞=%.2e (%s)\n", r.normF, r.note)
                    n_fail += 1
                    @printf(io, "%.2f,%.0f,%.4f,%.1f,%.6f,%s,%s,%.2f,%s\n",
                            gap_mm, Q_total, eps_org, r.delta_m * 1e6,
                            V_WARM, basename(r.path), false, r.wall_s, r.note)
                end
            catch err
                @printf("EXC: %s\n", err)
                n_fail += 1
                @printf(io, "%.2f,%.0f,%.4f,,,,%s,%.2f,exception:%s\n",
                        gap_mm, Q_total, eps_org, false, 0.0,
                        replace(string(err), "," => ";"))
            end
        end

        println()
        @printf("Built : %3d   Skipped: %3d   Failed: %3d\n",
                n_built, n_skip, n_fail)
        @printf("Index : %s\n", index_path)
    end
    println("=" ^ 72)
end

main()
