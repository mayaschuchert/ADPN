module Diffusivity

using ..Params

export D_mix, set_D_formulation!, D_formulation

# Module-level switch — lets the same residual chain use either the
# arithmetic mean (default) or the m_i-corrected D_eff (guide §4.2).
# Only changes behaviour in the two-phase regime (ε_org ≥ EPS_ORG_SAT);
# single-phase always returns D_aq regardless.
const _D_FORMULATION = Ref{Symbol}(:arithmetic)

"""
    set_D_formulation!(:arithmetic) or set_D_formulation!(:m_corrected)

`:arithmetic`   → D_mix = ε·D_org + (1−ε)·D_aq       (no m_i factor)
`:m_corrected`  → D_eff = ε·m_i·D_org + (1−ε)·D_aq   (partition-weighted)

For ions (D_org = 0) both give the same (1−ε)·D_aq.
For neutrals (AN, ADPN, PN) m_corrected gives much higher D at high ε.
"""
function set_D_formulation!(s::Symbol)
    @assert s in (:arithmetic, :m_corrected) "Unknown D formulation: $s"
    _D_FORMULATION[] = s
    return s
end

D_formulation() = _D_FORMULATION[]

"""
    D_mix(i, eps_org) -> Float64

Regime-aware effective diffusivity for species i.

**Single-phase regime** (ε_org < EPS_ORG_SAT ≈ 0.0794):
there are no organic droplets — all AN is dissolved in the aqueous phase.
Transport of every species happens through a single aqueous phase:
```
    D_i,mix = D_i,aq
```

**Two-phase regime** (ε_org ≥ EPS_ORG_SAT):
aqueous phase is saturated with AN; the remaining organic volume forms
droplets in local chemical equilibrium with the aqueous phase. The effective
volume-averaged diffusivity is the arithmetic mean (§4.1):
```
    D_i,mix = ε_org · D_i,org + (1 − ε_org) · D_i,aq
```
Ions have D_i,org = 0 (insoluble in organic), so for ions the two-phase
formula reduces to (1 − ε_org) · D_i,aq.

The guide's `D_eff` upgrade (§4.2) with partition coefficient m_i is a
drop-in replacement: `D_i,eff = ε_org · m_i · D_i,org + (1 − ε_org) · D_i,aq`
(same single-phase fallback).
"""
@inline function D_mix(i::Int, eps_org::Float64)
    if eps_org < EPS_ORG_SAT
        return D_aq[i]     # single-phase: no organic phase exists
    end
    # Two-phase regime — dispatch on the module-level formulation switch.
    mode = _D_FORMULATION[]
    if mode === :m_corrected
        return eps_org * m_partition[i] * D_org[i] + (1.0 - eps_org) * D_aq[i]
    else    # :arithmetic  (default)
        return eps_org * D_org[i] + (1.0 - eps_org) * D_aq[i]
    end
end

end # module
