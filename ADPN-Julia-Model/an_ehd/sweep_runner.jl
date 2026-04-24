# -----------------------------------------------------------------------------
# sweep_runner.jl — reusable pipeline for a single (ε_org, δ) sweep:
#   1. mesh + bulk equilibrium
#   2. initial guess, bootstrap (buffer ramp, kinetics ramp)
#   3. V continuation (AD Jacobian) + optional log-j continuation
#   4. cache per-V solutions, emit records/profile/meta CSVs
#
# Called by run_stage1.jl (single sweep) and run_stage2.jl (ε_org loop).
# -----------------------------------------------------------------------------
using Printf
using Dates

include(joinpath(@__DIR__, "ADPN_EHD.jl"))
using .ADPN_EHD
using .ADPN_EHD.Params
using .ADPN_EHD.Diffusivity
using .ADPN_EHD.Chemistry
using .ADPN_EHD.Kinetics
using .ADPN_EHD.Transport
using .ADPN_EHD.Assembly
using .ADPN_EHD.Solver
using .ADPN_EHD.Mesh

# Fixed, pipeline-wide constants
const N_MESH       = 100
const STRETCH      = 10.0
const V_START      = -1.0
const V_END        = -2.5
const V_BOOT       = V_START

const CACHE_DIR    = joinpath(@__DIR__, "output", "cache")
const DATA_DIR     = joinpath(@__DIR__, "output", "data")
isdir(CACHE_DIR) || mkpath(CACHE_DIR)
isdir(DATA_DIR)  || mkpath(DATA_DIR)

# -----------------------------------------------------------------------------
# Residual closure (captures mesh, ε_org, α's, c_eq)
# -----------------------------------------------------------------------------
residual_closure(mesh, eps_org, V, alpha_buf, alpha_kin, c_eq) =
    (F, u) -> full_residual!(F, u, mesh, eps_org, V, alpha_buf, alpha_kin, c_eq)

build_full_residual(mesh, eps_org, c_eq) =
    V -> residual_closure(mesh, eps_org, V, 1.0, 1.0, c_eq)

# -----------------------------------------------------------------------------
# Bootstrap — buffer ramp then kinetics ramp, parameterised by ε_org
# -----------------------------------------------------------------------------
function bootstrap!(u::Vector{Float64}, mesh, eps_org::Float64, c_eq;
                    verbose::Bool = true)

    # Step 1: zero-residual initial state (α_buf = α_kin = 0)
    residual! = residual_closure(mesh, eps_org, V_BOOT, 0.0, 0.0, c_eq)
    F0 = zeros(length(u))
    residual!(F0, u)
    normF = maximum(abs.(F0))
    verbose && @printf("[boot] α_buf=0, α_kin=0 at V=%.3f V:  |F|∞ = %.3e\n",
                       V_BOOT, normF)
    @assert normF < 1e-8 "Initial residual should be ~0 but got $(normF)"

    # Step 2: buffer ramp (α_buf 0 → 1 in 10 uniform steps)
    verbose && println("[boot] ramping α_buf 0 → 1 (10 uniform steps), α_kin = 0")
    for α in range(0.1, 1.0; length = 10)
        residual! = residual_closure(mesh, eps_org, V_BOOT, α, 0.0, c_eq)
        res = newton_solve!(u, residual!; max_iter = 25, verbose = false)
        verbose && @printf("  α_buf = %.2f   iter = %2d   |F|∞ = %.3e   %s\n",
                           α, res.iter, res.normF, res.converged ? "ok" : "FAIL")
        @assert res.converged "buffer ramp failed at α_buf = $α"
    end

    # Step 3: kinetics ramp (α_kin geometric 1e-6 → 1.0, 21 steps)
    verbose && println("[boot] ramping α_kin 1e-6 → 1.0 (geometric ×2), α_buf = 1")
    α_list = [min(1.0e-6 * 2.0^k, 1.0) for k in 0:20]
    α_list[end] = 1.0
    for α in α_list
        residual! = residual_closure(mesh, eps_org, V_BOOT, 1.0, α, c_eq)
        res = newton_solve!(u, residual!; max_iter = 25, verbose = false)
        verbose && @printf("  α_kin = %.3e  iter = %2d   |F|∞ = %.3e   %s\n",
                           α, res.iter, res.normF, res.converged ? "ok" : "FAIL")
        @assert res.converged "kinetics ramp failed at α_kin = $α"
    end

    verbose && println("[boot] completed at α_buf=1, α_kin=1, V = $(V_BOOT)")
    return u
end

# -----------------------------------------------------------------------------
# Physicality report — prints & returns a NamedTuple of metrics
# -----------------------------------------------------------------------------
function physicality_report(u::Vector{Float64}, mesh,
                            eps_org::Float64, delta_um::Float64,
                            V::Float64, c_eq;
                            alpha_buf = 1.0, alpha_kin = 1.0,
                            tag::String = "",
                            verbose::Bool = true)
    N = length(mesh.dx)
    c, phi = decode_state(u, N)

    pH_bulk = -log10(c[1, N] / 1000.0)
    phi_span = maximum(phi) - minimum(phi)

    eneutr_9sp = zeros(N)
    Na_dev     = zeros(N)
    for ix in 1:N
        c_Na_loc = c[2, ix] + c[3, ix] + 2*c[4, ix] + 3*c[5, ix] - c[1, ix]
        s9       = c[1,ix] - c[2,ix] - c[3,ix] - 2*c[4,ix] - 3*c[5,ix] + c_Na_loc
        eneutr_9sp[ix] = abs(s9)
        Na_dev[ix]     = abs(c_Na_loc - c_eq.Na) / c_eq.Na
    end

    J_charge_face = zeros(N + 1)
    for ix in 1:(N - 1)
        dx_face = 0.5 * (mesh.dx[ix] + mesh.dx[ix + 1])
        s = 0.0
        for k in 1:5
            D_k = D_mix(k, eps_org)
            s += z_species[k] * sg_flux(c[k, ix], c[k, ix + 1],
                                        phi[ix], phi[ix + 1],
                                        D_k, z_species[k], dx_face)
        end
        J_charge_face[ix + 1] = s
    end
    j1, j2, j3 = tafel_currents(c[6, 1], phi[1], V, alpha_kin)
    J_charge_face[1] = -(j1 + j2 + j3) / Params.F
    iface_var = maximum(abs.(diff(J_charge_face[2:N])))

    R_buf = zeros(8)
    buffer_sources!(R_buf, c[1,1], c[2,1], c[3,1], c[4,1], c[5,1], alpha_buf)
    max_Rbuf_surf = maximum(abs.(R_buf[1:5]))
    buffer_sources!(R_buf, c[1,N], c[2,N], c[3,N], c[4,N], c[5,N], alpha_buf)
    max_Rbuf_bulk = maximum(abs.(R_buf[1:5]))

    c_AN_bulk_val = c[6, N]
    c_AN_surface  = c[6, 1]
    AN_depletion  = c_AN_surface / c_AN_bulk_val

    j_tot   = j1 + j2 + j3
    FE_ADPN = 100.0 * j1 / (j_tot + 1e-300)
    FE_PN   = 100.0 * j2 / (j_tot + 1e-300)
    FE_HER  = 100.0 * j3 / (j_tot + 1e-300)

    if verbose
        println("────────── PHYSICALITY  ", tag, " ──────────")
        @printf("  V = %+7.3f V vs SHE,  ε_org = %.3f,  δ = %.1f μm\n",
                V, eps_org, delta_um)
        @printf("  j₁=%.3e  j₂=%.3e  j₃=%.3e  [A/m²]\n", j1, j2, j3)
        @printf("  j_total = %.3e A/m² = %.2f mA/cm²\n", j_tot, j_tot * 0.1)
        @printf("  FE_ADPN = %6.2f %%   FE_PN = %6.2f %%   FE_HER = %6.2f %%\n",
                FE_ADPN, FE_PN, FE_HER)
        @printf("  bulk pH                   : %.3f\n", pH_bulk)
        @printf("  φ_l span across DL        : %.3e V\n", phi_span)
        @printf("  max |Σ z_i c_i| (9-sp)    : %.3e mol m⁻³\n", maximum(eneutr_9sp))
        @printf("  max Na⁺ DL dev            : %.3e\n", maximum(Na_dev))
        @printf("  Σ z_i N_i interior  var   : %.3e mol m⁻² s⁻¹\n", iface_var)
        @printf("  max |R_buf| @ surface     : %.3e mol m⁻³ s⁻¹\n", max_Rbuf_surf)
        @printf("  max |R_buf| @ bulk        : %.3e mol m⁻³ s⁻¹\n", max_Rbuf_bulk)
        @printf("  c_AN(0) / c_AN,bulk       : %.4f\n", AN_depletion)
        println("──────────────────────────────────────────────")
    end

    return (
        V = V, eps_org = eps_org, delta_um = delta_um,
        j1 = j1, j2 = j2, j3 = j3, j_total = j_tot,
        FE_ADPN = FE_ADPN, FE_PN = FE_PN, FE_HER = FE_HER,
        pH_bulk = pH_bulk, phi_span = phi_span,
        max_eneutr = maximum(eneutr_9sp), Na_dev = maximum(Na_dev),
        iface_var = iface_var,
        max_Rbuf_surf = max_Rbuf_surf, max_Rbuf_bulk = max_Rbuf_bulk,
        AN_depletion = AN_depletion,
        c_AN_surface = c_AN_surface, c_AN_bulk = c_AN_bulk_val,
        phi_l_surface = phi[1],
    )
end

# -----------------------------------------------------------------------------
# Caching & CSV export
# -----------------------------------------------------------------------------
function save_solution(u::Vector{Float64}, eps_org::Float64, delta_um::Real, V::Float64)
    fn = @sprintf("s_eo%.3f_d%d_V%.6f.bin", eps_org, round(Int, delta_um), V)
    path = joinpath(CACHE_DIR, fn)
    open(path, "w") do f
        write(f, Int64(length(u)))
        write(f, u)
    end
    return path
end

function export_records_csv(records, stage_tag::String, eps_org, delta_um)
    path = joinpath(DATA_DIR,
        @sprintf("%s_records_eo%.3f_d%d.csv",
                 stage_tag, eps_org, round(Int, delta_um)))
    open(path, "w") do io
        println(io, "V,eps_org,delta_um,j1,j2,j3,j_total,FE_ADPN,FE_PN,FE_HER,",
                     "pH_bulk,phi_span,max_eneutr,Na_dev,iface_var,max_Rbuf_surf,",
                     "max_Rbuf_bulk,AN_depletion,c_AN_surface,c_AN_bulk,phi_l_surface")
        for r in records
            @printf(io, "%.8e,%.4f,%.4f,%.8e,%.8e,%.8e,%.8e,%.6f,%.6f,%.6f,",
                    r.V, r.eps_org, r.delta_um, r.j1, r.j2, r.j3, r.j_total,
                    r.FE_ADPN, r.FE_PN, r.FE_HER)
            @printf(io, "%.6f,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.6f,%.8e,%.8e,%.8e\n",
                    r.pH_bulk, r.phi_span, r.max_eneutr, r.Na_dev, r.iface_var,
                    r.max_Rbuf_surf, r.max_Rbuf_bulk, r.AN_depletion,
                    r.c_AN_surface, r.c_AN_bulk, r.phi_l_surface)
        end
    end
    println("[export] records → $path")
end

function export_profile_at(history, mesh, stage_tag::String,
                           eps_org, delta_um, V_target::Real)
    V_target_f = Float64(V_target)
    best_i = argmin([abs(h[1] - V_target_f) for h in history])
    V, u = history[best_i]
    N = length(mesh.dx)
    c, phi = decode_state(u, N)
    path = joinpath(DATA_DIR,
        @sprintf("%s_profile_eo%.3f_d%d_V%.3f.csv",
                 stage_tag, eps_org, round(Int, delta_um), V))
    open(path, "w") do io
        println(io, "x_m,c_H,c_OH,c_H2PO4,c_HPO4,c_PO4,c_AN,c_ADPN,c_PN,phi_l")
        for ix in 1:N
            @printf(io, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                    mesh.x_c[ix], c[1,ix], c[2,ix], c[3,ix], c[4,ix],
                    c[5,ix], c[6,ix], c[7,ix], c[8,ix], phi[ix])
        end
    end
    println("[export] profile @ V=$(V) (target $V_target_f) → $path")
end

function export_meta(c_eq, mesh, stage_tag::String, eps_org, delta_um)
    path = joinpath(DATA_DIR,
        @sprintf("%s_meta_eo%.3f_d%d.txt", stage_tag, eps_org, round(Int, delta_um)))
    open(path, "w") do io
        println(io, "stage=", stage_tag)
        println(io, "eps_org=", eps_org)
        println(io, "delta_um=", delta_um)
        println(io, "N_mesh=", length(mesh.dx))
        println(io, "stretch=", STRETCH)
        println(io, "pH_bulk=", c_eq.pH)
        println(io, "cH_bulk=", c_eq.H)
        println(io, "cOH_bulk=", c_eq.OH)
        println(io, "cH2PO4_bulk=", c_eq.H2PO4)
        println(io, "cHPO4_bulk=", c_eq.HPO4)
        println(io, "cPO4_bulk=", c_eq.PO4)
        println(io, "cNa_bulk=", c_eq.Na)
        println(io, "cAN_bulk=", c_AN_bulk(eps_org))
        println(io, "dx_min=", minimum(mesh.dx))
        println(io, "dx_max=", maximum(mesh.dx))
    end
    println("[export] meta → $path")
end

# -----------------------------------------------------------------------------
# run_one_sweep — the full pipeline for a single (ε_org, δ)
# -----------------------------------------------------------------------------
function run_one_sweep(eps_org::Float64, delta_um::Float64;
                       stage_tag::String = "stage1",
                       profile_V_targets = (-1.5, -1.8, -2.0, -2.5),
                       verbose::Bool = true)

    delta = delta_um * 1.0e-6

    verbose && println("════════════════════════════════════════════════════════════")
    verbose && @printf(" %s — ε_org = %.3f, δ = %.1f μm\n", stage_tag, eps_org, delta_um)
    verbose && println(" Start: ", Dates.now())
    verbose && println("════════════════════════════════════════════════════════════")

    mesh = make_mesh(N_MESH, delta; stretch = STRETCH)
    verbose && @printf("[mesh]  N=%d dx_min=%.3e dx_max=%.3e ratio=%.2f\n",
                       N_MESH, minimum(mesh.dx), maximum(mesh.dx),
                       maximum(mesh.dx)/minimum(mesh.dx))

    c_eq = solve_phosphate_equilibrium()
    verbose && @printf("[eq]   pH=%.3f cOH=%.3e cHPO4=%.3e cPO4=%.3e cNa=%.1f  cAN,bulk=%.1f\n",
                       c_eq.pH, c_eq.OH, c_eq.HPO4, c_eq.PO4, c_eq.Na, c_AN_bulk(eps_org))

    # Equilibrium sanity
    @assert abs(c_eq.H * c_eq.OH - K_w) / K_w < 1e-8
    @assert abs(c_eq.H * c_eq.HPO4 / c_eq.H2PO4 - K_a2) / K_a2 < 1e-8
    @assert abs(c_eq.H * c_eq.PO4  / c_eq.HPO4  - K_a3) / K_a3 < 1e-8
    R_chk = zeros(8)
    buffer_sources!(R_chk, c_eq.H, c_eq.OH, c_eq.H2PO4, c_eq.HPO4, c_eq.PO4)
    @assert maximum(abs.(R_chk[1:5])) < 1e-6

    # Initial guess & bootstrap
    u = make_initial_guess(N_MESH, c_eq, eps_org)
    bootstrap!(u, mesh, eps_org, c_eq; verbose = verbose)

    if verbose
        physicality_report(u, mesh, eps_org, delta_um, V_START, c_eq;
                           tag = "post-bootstrap")
    end

    # V continuation (AD Jacobian — most robust near rank-deficiency)
    verbose && println("\n[cont] Phase A: V continuation from $(V_START) V (AD Jacobian)")
    hist_V = newton_continuation(u, V_START, V_END,
                                 build_full_residual(mesh, eps_org, c_eq);
                                 ds_init = 0.05, ds_min = 1e-4, ds_max = 0.20,
                                 max_iter = 60,
                                 jacobian_mode = :ad,
                                 verbose = false)
    V_last = hist_V[end][1]
    verbose && @printf("[cont] Phase A stopped at V = %.4f V  (%d points)\n",
                       V_last, length(hist_V))

    # log-j fallback if V continuation stalls before V_END
    history = copy(hist_V)
    if V_last > V_END + 1e-6
        idx_restart = argmin([abs(h[1] - (-2.0)) for h in hist_V])
        V_restart, u_restart = hist_V[idx_restart]
        verbose && @printf("[cont] Phase B: log-j restart at V=%.4f (idx %d)\n",
                           V_restart, idx_restart)
        faradaic_current = (u_state, V_now) -> begin
            c_AN_s = exp(clamp(u_state[6], -50.0, 50.0))
            phi_s  = u_state[9]
            j1, j2, j3 = tafel_currents(c_AN_s, phi_s, V_now, 1.0)
            return j1 + j2 + j3
        end
        hist_logj = newton_continuation_logj(u_restart, V_restart, V_END,
                                             build_full_residual(mesh, eps_org, c_eq),
                                             faradaic_current;
                                             n_steps_per_decade = 30,
                                             max_iter = 40, verbose = false)
        history = vcat(hist_V[1:idx_restart], hist_logj[2:end])
        verbose && @printf("[cont] Phase B reached V=%.4f (+%d log-j points)\n",
                           hist_logj[end][1], length(hist_logj) - 1)
    end

    verbose && @printf("[cont] recorded %d converged points, V∈[%.3f, %.3f]\n",
                       length(history),
                       minimum(first.(history)), maximum(first.(history)))

    # Cache + per-point physicality
    all_records = NamedTuple[]
    for (V, u_V) in history
        save_solution(u_V, eps_org, delta_um, V)
        rec = physicality_report(u_V, mesh, eps_org, delta_um, V, c_eq;
                                 tag = @sprintf("V=%.3f", V),
                                 verbose = false)
        push!(all_records, rec)
    end

    # Export
    export_records_csv(all_records, stage_tag, eps_org, delta_um)
    for V_target in profile_V_targets
        export_profile_at(history, mesh, stage_tag, eps_org, delta_um, V_target)
    end
    export_meta(c_eq, mesh, stage_tag, eps_org, delta_um)

    verbose && println("\n════════════════════════════════════════════════════════════")
    verbose && @printf(" %s for ε_org=%.3f complete.  End: %s\n",
                       stage_tag, eps_org, Dates.now())
    verbose && println("════════════════════════════════════════════════════════════")

    return history, all_records
end
