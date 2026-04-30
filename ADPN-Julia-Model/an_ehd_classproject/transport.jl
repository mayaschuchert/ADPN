module Transport

using ..Params

export sg_flux

"""
    sg_flux(c_L, c_R, phi_L, phi_R, D_mix, z_i, dx) -> Float64

Scharfetter–Gummel face flux N_i [mol m⁻² s⁻¹], oriented +x (from L to R).

  α   = z_i F (φ_R − φ_L) / (R T)
  N_i = − D_mix/dx · [ B(α) c_R − B(−α) c_L ]
  B(α) = α / (e^α − 1)

For z = 0 (or |α| < 1e−10) reduces to − D_mix (c_R − c_L) / dx.
α is clamped to [−700, 700] to avoid float overflow.
"""
@inline function sg_flux(c_L::Real, c_R::Real,
                         phi_L::Real, phi_R::Real,
                         D::Real, z_i::Int, dx::Real)

    if z_i == 0
        return -D * (c_R - c_L) / dx
    end

    alpha = z_i * F * (phi_R - phi_L) / (R_gas * T)
    alpha = clamp(alpha, -700.0, 700.0)

    # Bernoulli functions B(α) = α / (e^α − 1), B(−α) = α·e^α / (e^α − 1).
    # The naïve form divides by zero as α → 0. Use a Taylor expansion for
    # small |α| so the function remains smooth AND its derivative w.r.t. φ
    # is continuous through α = 0. This is essential for ForwardDiff AD
    # (a hard if/else branch with a centered-difference fallback gives a
    # spurious zero derivative exactly at the α = 0 seed point).
    if abs(alpha) < 0.01
        # Taylor series: B(α)  = 1 − α/2 + α²/12 − α⁴/720 + ...
        #                 B(−α) = 1 + α/2 + α²/12 + α⁴/720 + ...
        a2 = alpha * alpha
        B_pos = 1.0 - alpha / 2.0 + a2 / 12.0
        B_neg = 1.0 + alpha / 2.0 + a2 / 12.0
    else
        ea  = exp(alpha)
        den = ea - 1.0
        B_pos = alpha / den
        B_neg = alpha * ea / den
    end
    return -D / dx * (B_pos * c_R - B_neg * c_L)
end

end # module
