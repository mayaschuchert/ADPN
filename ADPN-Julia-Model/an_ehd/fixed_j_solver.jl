module FixedJ

# v6 §20 — bisect V vs SHE so the model's total Faradaic current density
# matches a target j (mA/cm² → A/m² conversion happens at the call site).
# Inside, this wraps the existing v5 Newton continuation: at each candidate
# V the residual is built and Newton is run from a warm start.

using ..Params
using ..Mesh
using ..Chemistry
using ..Kinetics
using ..Assembly
using ..Solver

export solve_at_j, FixedJResult

# ---------- Result type ----------
struct FixedJResult
    converged::Bool
    V_cathode::Float64        # V vs SHE
    j_total::Float64          # A m⁻², achieved
    j1::Float64               # ADPN partial current density
    j2::Float64               # PN partial current density
    j3::Float64               # HER partial current density
    FE_ADN_pct::Float64
    FE_PN_pct::Float64
    FE_HER_pct::Float64
    state::Vector{Float64}    # converged DOF vector
    n_bisect::Int
    note::String
end

# ---------- Newton solve at one V using a warm start ----------
function _solve_at_V!(u::Vector{Float64}, mesh, eps_org::Float64,
                      V::Float64, c_eq;
                      max_iter::Int = 60,
                      jacobian_mode::Symbol = :ad,
                      verbose::Bool = false)
    residual! = (F, x) -> full_residual!(F, x, mesh, eps_org, V,
                                         1.0, 1.0, c_eq)
    return newton_solve!(u, residual!;
                         max_iter      = max_iter,
                         jacobian_mode = jacobian_mode,
                         verbose       = verbose)
end

# ---------- Total Faradaic current density at the converged state ----------
function _j_total(u::Vector{Float64}, V::Float64)
    j1, j2, j3 = faradaic_currents_from_state(u, V, 1.0)
    return (j1 + j2 + j3, j1, j2, j3)
end

"""
    solve_at_j(j_target_A_m2, eps_org, delta_m, mesh, u_warm, c_eq;
               V_lo=-2.5, V_hi=-0.8, tol_rel=1e-3,
               max_bisect=40, jacobian_mode=:ad,
               j0=nothing, alpha_c=nothing, verbose=false) -> FixedJResult

Bisect V vs SHE in [V_lo, V_hi] (negative range; V_lo more negative) until the
model's total cathodic current density matches `j_target_A_m2` to within
`tol_rel · j_target`.

`u_warm` is the warm-start DOF vector (must be a converged solution at any V
in the bracket — typically pulled from the Stage 3 cache for the matching
(eps_org, delta) tuple). The function does **not** mutate `u_warm`; it copies
internally.

If `j0` and `alpha_c` are provided as 3-tuples, they are pushed into the
Kinetics override Ref for the duration of this call (restored on exit). This
is how Stage 4 evaluates the loss at trial parameter sets without touching
Params constants. Pass `nothing` to use Params defaults (also default).
"""
function solve_at_j(j_target_A_m2::Float64,
                    eps_org::Float64,
                    delta_m::Float64,
                    mesh,
                    u_warm::Vector{Float64},
                    c_eq;
                    V_lo::Float64       = -2.5,
                    V_hi::Float64       = -0.8,
                    tol_rel::Float64    = 1e-3,
                    max_bisect::Int     = 40,
                    jacobian_mode::Symbol = :ad,
                    newton_max_iter::Int  = 60,
                    j0::Union{Nothing,NTuple{3,Float64}}      = nothing,
                    alpha_c::Union{Nothing,NTuple{3,Float64}} = nothing,
                    verbose::Bool       = false)

    @assert V_lo < V_hi "V_lo must be more negative than V_hi"
    @assert j_target_A_m2 > 0 "j_target must be positive (cathodic)"

    # Optional kinetics override (restored on exit)
    use_override = !(j0 === nothing || alpha_c === nothing)
    prev_override = Kinetics.KIN_OVERRIDE[]
    if use_override
        set_kinetic_override!(j0, alpha_c)
    end

    try
        # Bracket evaluation: solve at V_hi (less negative, low j) and V_lo
        u_hi = copy(u_warm)
        res_hi = _solve_at_V!(u_hi, mesh, eps_org, V_hi, c_eq;
                              max_iter = newton_max_iter,
                              jacobian_mode = jacobian_mode,
                              verbose = verbose)
        if !res_hi.converged
            return FixedJResult(false, V_hi, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                                u_warm, 0, "Newton failed at V_hi=$(V_hi)")
        end
        j_hi, _, _, _ = _j_total(u_hi, V_hi)

        u_lo = copy(u_warm)
        res_lo = _solve_at_V!(u_lo, mesh, eps_org, V_lo, c_eq;
                              max_iter = newton_max_iter,
                              jacobian_mode = jacobian_mode,
                              verbose = verbose)
        if !res_lo.converged
            return FixedJResult(false, V_lo, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                                u_warm, 0, "Newton failed at V_lo=$(V_lo)")
        end
        j_lo, _, _, _ = _j_total(u_lo, V_lo)

        if !(j_hi < j_target_A_m2 < j_lo)
            return FixedJResult(false, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                                u_warm, 0,
                                "j_target=$(j_target_A_m2) outside bracket " *
                                "[j(V_hi)=$(j_hi), j(V_lo)=$(j_lo)]")
        end

        # Bisection on V (monotone j(V) in this range — verified at startup).
        # Note: V_lo is more negative (higher |V|, higher j), V_hi less negative.
        # Use the converged states at the bracket as warm starts to keep
        # Newton in-basin throughout the bisection.
        a, b = V_lo, V_hi          # j(a) > j_target > j(b)
        u_a, u_b = u_lo, u_hi
        j_a, j_b = j_lo, j_hi

        u_mid = copy(u_b)
        V_mid = 0.5 * (a + b)
        j_mid = NaN
        j1m = j2m = j3m = NaN

        n_bisect = 0
        for k in 1:max_bisect
            V_mid = 0.5 * (a + b)
            # Pick the closer warm start to V_mid as initial guess
            u_seed = abs(V_mid - a) < abs(V_mid - b) ? u_a : u_b
            u_mid  = copy(u_seed)
            res_m  = _solve_at_V!(u_mid, mesh, eps_org, V_mid, c_eq;
                                  max_iter = newton_max_iter,
                                  jacobian_mode = jacobian_mode,
                                  verbose = verbose)
            if !res_m.converged
                return FixedJResult(false, V_mid, NaN, NaN, NaN, NaN, NaN, NaN, NaN,
                                    u_warm, k,
                                    "Newton failed at V_mid=$(V_mid) (bisect step $k)")
            end
            j_mid, j1m, j2m, j3m = _j_total(u_mid, V_mid)
            if abs(j_mid - j_target_A_m2) <= tol_rel * j_target_A_m2
                n_bisect = k
                FE_total = j1m + j2m + j3m
                FE_ADN  = 100 * j1m / FE_total
                FE_PN   = 100 * j2m / FE_total
                FE_HER  = 100 * j3m / FE_total
                return FixedJResult(true, V_mid, j_mid, j1m, j2m, j3m,
                                    FE_ADN, FE_PN, FE_HER,
                                    u_mid, n_bisect, "")
            end
            if j_mid > j_target_A_m2
                a = V_mid; u_a = u_mid; j_a = j_mid
            else
                b = V_mid; u_b = u_mid; j_b = j_mid
            end
            n_bisect = k
        end

        # Did not converge to tol_rel within max_bisect; return best-effort
        FE_total = j1m + j2m + j3m
        FE_ADN  = 100 * j1m / FE_total
        FE_PN   = 100 * j2m / FE_total
        FE_HER  = 100 * j3m / FE_total
        return FixedJResult(false, V_mid, j_mid, j1m, j2m, j3m,
                            FE_ADN, FE_PN, FE_HER,
                            u_mid, n_bisect,
                            "Bisection did not reach tol_rel=$(tol_rel) in " *
                            "$(max_bisect) steps; |Δj/j| ≈ " *
                            "$(abs(j_mid - j_target_A_m2)/j_target_A_m2)")
    finally
        if use_override
            Kinetics.KIN_OVERRIDE[] = prev_override
        end
    end
end

end # module
