module Chemistry

using ..Params

# Lightweight bisection (avoids dependency on Roots.jl)
function _bisect(f::F, a::Float64, b::Float64;
                 tol::Float64 = 1e-14, max_iter::Int = 200) where {F<:Function}
    fa = f(a); fb = f(b)
    if fa == 0.0;  return a; end
    if fb == 0.0;  return b; end
    @assert sign(fa) != sign(fb) "bisection: sign(f(a)) == sign(f(b)) in [$a, $b]"
    for _ in 1:max_iter
        m  = 0.5 * (a + b)
        fm = f(m)
        if fm == 0.0 || 0.5 * (b - a) < tol
            return m
        end
        if sign(fm) == sign(fa)
            a = m; fa = fm
        else
            b = m; fb = fm
        end
    end
    return 0.5 * (a + b)
end

export solve_phosphate_equilibrium,
       buffer_sources!,
       make_initial_guess,
       bulk_concentration,
       c_AN_bulk

# ------------------------------------------------------------
# 6.4  Bulk equilibrium (phosphate charge balance)
# ------------------------------------------------------------
"""
    solve_phosphate_equilibrium(C_P_total, c_Na_input) -> NamedTuple

Solve bulk phosphate + water charge balance in log10(c_H). All conc. in mol m⁻³.
Returns (H, OH, H2PO4, HPO4, PO4, Na, pH).
"""
function solve_phosphate_equilibrium(C_P_total::Float64 = Params.C_P_total,
                                     c_Na_input::Float64 = Params.c_Na_input)
    function charge_balance(log10_cH)
        cH     = 10.0^log10_cH
        cOH    = K_w / cH
        denom  = cH / K_a2 + 1.0 + K_a3 / cH
        cHPO4  = C_P_total / denom
        cH2PO4 = cH * cHPO4 / K_a2
        cPO4   = K_a3 * cHPO4 / cH
        return c_Na_input + cH - (cOH + cH2PO4 + 2.0*cHPO4 + 3.0*cPO4)
    end

    # log10(cH [mol/m³]) brackets: pH 6..15 in std units → cH = 10⁻³ .. 10⁻¹² mol/m³
    log10_cH_eq = _bisect(charge_balance, -12.0, -3.0)

    cH     = 10.0^log10_cH_eq
    cOH    = K_w / cH
    denom  = cH / K_a2 + 1.0 + K_a3 / cH
    cHPO4  = C_P_total / denom
    cH2PO4 = cH * cHPO4 / K_a2
    cPO4   = K_a3 * cHPO4 / cH
    pH_std = -log10(cH / 1000.0)

    return (H = cH, OH = cOH, H2PO4 = cH2PO4, HPO4 = cHPO4, PO4 = cPO4,
            Na = c_Na_input, pH = pH_std)
end

# ------------------------------------------------------------
# 6.2  Buffer source terms R_buf,i
# ------------------------------------------------------------
"""
    buffer_sources!(R, c_H, c_OH, c_H2PO4, c_HPO4, c_PO4, alpha_buf=1.0)

Fill R::AbstractVector (length ≥ 8) with buffer rates [mol m⁻³ s⁻¹]:
indices 1..5 are H⁺, OH⁻, H₂PO₄⁻, HPO₄²⁻, PO₄³⁻; 6..8 (AN/ADPN/PN) set to 0.

**OH⁻-pathway formulation**:
  R1:  H₂O  ⇌  H⁺ + OH⁻
  R2:  OH⁻ + H₂PO₄⁻  ⇌  HPO₄²⁻ + H₂O
  R3:  OH⁻ + HPO₄²⁻  ⇌  PO₄³⁻  + H₂O

Net sources per species (stoichiometry table in guide discussion):
  S_H      =  r₁
  S_OH     =  r₁ − r₂ − r₃
  S_H₂PO₄  = −r₂
  S_HPO₄   =  r₂ − r₃
  S_PO₄    =  r₃

Σ zᵢ·Sᵢ = 0 for charge conservation. Equilibrium fixed-point coincides with
the H⁺-pathway form but kinetic timescales differ — OH⁻ is now consumed
directly by r₂ and r₃ at high pH rather than via slow water autoprotolysis.
"""
function buffer_sources!(R::AbstractVector,
                         c_H, c_OH, c_H2PO4, c_HPO4, c_PO4,
                         alpha_buf::Float64 = 1.0)
    # R1 — water autoprotolysis (unchanged form, updated k₁f/k₁r)
    r1 = k1f - k1r * c_H * c_OH
    # R2 — OH⁻ + H₂PO₄⁻  ⇌  HPO₄²⁻ + H₂O
    r2 = k2f * c_OH * c_H2PO4 - k2r * c_HPO4
    # R3 — OH⁻ + HPO₄²⁻  ⇌  PO₄³⁻  + H₂O
    r3 = k3f * c_OH * c_HPO4  - k3r * c_PO4

    R[1] = alpha_buf * ( r1              )   # H⁺     — only from water
    R[2] = alpha_buf * ( r1 - r2 - r3    )   # OH⁻    — water + consumed by r₂, r₃
    R[3] = alpha_buf * (-r2              )   # H₂PO₄⁻ — consumed by r₂
    R[4] = alpha_buf * ( r2 - r3         )   # HPO₄²⁻ — produced by r₂, consumed by r₃
    R[5] = alpha_buf * ( r3              )   # PO₄³⁻  — produced by r₃
    R[6] = 0.0                                # AN
    R[7] = 0.0                                # ADPN
    R[8] = 0.0                                # PN
    return R
end

# ------------------------------------------------------------
# 7.2  Bulk AN regime
# ------------------------------------------------------------
"""
    c_AN_bulk(eps_org) -> Float64  [mol m⁻³]

Bulk aqueous AN concentration for a given loading `eps_org` (volume fraction
of AN added per total volume).

**Convention A — per total solution volume.**

**Single-phase regime** (ε_org < EPS_ORG_SAT ≈ 0.0862): all AN dissolves in
the aqueous phase. The homogeneous single-phase mixture has total volume
V_total (≈ V_water + V_AN when volumes add ideally), and the concentration
of dissolved AN is moles-of-AN divided by total volume:
```
c_AN,aq = ε_org · (ρ_AN / M_AN)
```
Linear in ε_org, zero at ε_org = 0 (no AN — pathological limit).

**Two-phase regime** (ε_org ≥ EPS_ORG_SAT): aqueous phase is saturated at
C_AN_SAT = ρ_AN/(M_AN · m_AN). Excess AN exists as organic droplets. The
aqueous concentration is pinned at the thermodynamic saturation value.

Continuous at ε_org = EPS_ORG_SAT where ε_sat · ρ_AN/M_AN = C_AN_SAT exactly.
"""
@inline function c_AN_bulk(eps_org::Float64; eps_sat::Float64 = EPS_ORG_SAT)
    if eps_org < eps_sat
        return eps_org * MOLAR_DENSITY_AN
    else
        return C_AN_SAT
    end
end

"""
    bulk_concentration(k, c_eq, eps_org) -> Float64

Dirichlet bulk conc. for species k. ADPN/PN seeded at 1e-3 mol m⁻³.
"""
@inline function bulk_concentration(k::Int, c_eq, eps_org::Float64)
    k == 1 && return c_eq.H
    k == 2 && return c_eq.OH
    k == 3 && return c_eq.H2PO4
    k == 4 && return c_eq.HPO4
    k == 5 && return c_eq.PO4
    k == 6 && return c_AN_bulk(eps_org)
    k == 7 && return 1.0e-3
    k == 8 && return 1.0e-3
    error("bulk_concentration: invalid species $k")
end

# ------------------------------------------------------------
# 6.4  Initial guess builder
# ------------------------------------------------------------
"""
    make_initial_guess(N_mesh, c_eq, eps_org; c_seed=1e-3) -> Vector{Float64}

Build 9·N vector of DOFs (cell-major): 8 log-concentrations + φ_l per cell.
With α_buf = α_kin = 0 the residual is identically zero for this guess.
"""
function make_initial_guess(N_mesh::Int, c_eq, eps_org::Float64;
                            c_seed::Float64 = 1e-3)
    u0 = zeros(9 * N_mesh)
    log_c0 = [log(c_eq.H),    log(c_eq.OH),
              log(c_eq.H2PO4), log(c_eq.HPO4), log(c_eq.PO4),
              log(c_AN_bulk(eps_org)),
              log(c_seed),     # ADPN seed
              log(c_seed)]     # PN seed
    @assert length(log_c0) == 8
    for ix in 1:N_mesh
        base = 9 * (ix - 1)
        @inbounds for k in 1:8
            u0[base + k] = log_c0[k]
        end
        u0[base + 9] = 0.0   # φ_l = 0 everywhere
    end
    return u0
end

end # module
