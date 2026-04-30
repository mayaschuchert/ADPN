module Mesh

export make_mesh

"""
    make_mesh(N, delta; stretch=10.0) -> NamedTuple

Cell-centred FV mesh on [0, δ] with geometric grading. N cells, finest at x = 0.
`stretch` = dx_max / dx_min. Returns (dx, x_c, x_faces).

  dx[k]   = dx_min · r^(k-1)           k = 1..N,   r = stretch^(1/(N-1))
  dx_min  = δ (r − 1) / (r^N − 1)

Fills [0, δ] exactly for any δ.
"""
function make_mesh(N::Int, delta::Float64; stretch::Float64 = 10.0)
    @assert N ≥ 2        "N must be ≥ 2"
    @assert delta > 0    "delta must be positive"
    @assert stretch ≥ 1  "stretch must be ≥ 1"

    r      = stretch^(1.0 / (N - 1))
    dx_min = delta * (r - 1.0) / (r^N - 1.0)
    dx     = [dx_min * r^(k - 1) for k in 1:N]

    # cell-centre positions:
    x_c = cumsum(dx) .- dx ./ 2.0

    # face positions (N+1):  face 1 = 0 (electrode), face N+1 = δ (bulk).
    x_faces = vcat(0.0, cumsum(dx))

    return (dx = dx, x_c = x_c, x_faces = x_faces)
end

end # module
