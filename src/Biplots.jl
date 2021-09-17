module Biplots

using LinearAlgebra
using Statistics

import Makie

"""
    biplot(X)

Biplot of design matrix `X`.

# Biplot attributes

* `dim` - number of dimensions `dim ∈ {2,3}` (default to `2`)
* `f`   - transformation function (default to `f(X) = X .- mean(X, dims=1)`)
* `α`   - shape parameter `α ∈ [0,1]` (default to `1`)
* `κ`   - normalization constant for axes (default to `√(size(X,1)-1)`)

# Aesthetics attributes

* `axesbody`  - size of principal axes' body
* `axeshead`  - size of principal axes' head
* `axescolor` - color of principal axes
* `axesnames` - names of principal axes
* `dotsize`   - size of samples
* `dotcolor`  - color of samples

See https://en.wikipedia.org/wiki/Biplot.

## References

* Gabriel, K. 1971. [The Biplot Graphic Display of Matrices
  with Application to Principal Component Analysis]
  (https://www.jstor.org/stable/2334381)

* Aitchison, J. & Greenacre, M. 2002. [Biplots of Compositional data]
  (https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9876.00275)
"""
@Makie.recipe(Biplot, X) do scene
  Makie.Attributes(;
    # biplot attributes
    dim = 2,
    f   = X -> X .- mean(X, dims=1),
    α   = 1.0,
    κ   = nothing,

    # aesthetic attributes
    axesbody  = nothing,
    axeshead  = nothing,
    axescolor = :black,
    axesnames = nothing,
    dotsize   = nothing,
    dotcolor  = :black,
    dotnames  = nothing,
    showdotnames = true,
  )
end

function Makie.plot!(plot::Biplot{<:Tuple{AbstractMatrix}})
  # retrieve parameters
  X = plot[:X][]
  d = plot[:dim][]
  f = plot[:f][]
  α = plot[:α][]
  κ = plot[:κ][]

  # size of design matrix
  n, m = size(X)

  # default normalization of axes
  if isnothing(κ)
    κ = α > 0.5 ? √(n-1) : 1.0
  end

  # retrieve options
  axesbody     = plot[:axesbody][]
  axeshead     = plot[:axeshead][]
  axescolor    = plot[:axescolor][]
  axesnames    = plot[:axesnames][]
  dotsize      = plot[:dotsize][]
  dotcolor     = plot[:dotcolor][]
  dotnames     = plot[:dotnames][]
  showdotnames = plot[:showdotnames][]

  # defaults differ on 2 or 3 dimensions
  if isnothing(axesbody)
    axesbody = d == 2 ? 2 : 0.01
  end
  if isnothing(axeshead)
    axeshead = d == 2 ? 6 : 0.03
  end
  if isnothing(axesnames)
    axesnames = ["x$i" for i in 1:m]
  end
  if isnothing(dotsize)
    dotsize = d == 2 ? 4 : 10
  end
  if isnothing(dotnames)
    dotnames = ["i" for i in 1:n]
  end

  # sanity checks
  @assert d ∈ [2,3] "d must be 2 or 3"
  @assert 0 ≤ α ≤ 1 "α must be in [0,1]"
  @assert length(axesnames) == m "axesnames must have length $m"
  @assert length(dotnames) == n "dotnames must have length $n"

  # transformation
  Z = f(X)

  # singular value decomposition
  U, σ, V = svd(Z)

  # variance explained
  v = σ[1:d] .^ 2 / sum(σ .^ 2)

  # matrix factors X ≈ F*G'
  F = U[:,1:d] .* (σ[1:d] .^ α)'
  G = V[:,1:d] .* (σ[1:d] .^ (1 - α))'

  # plot principal axes
  points = fill(Makie.Point(ntuple(i->0., d)), n)
  direcs = [Makie.Vec{d}(v) for v in eachrow(G)] ./ κ
  Makie.arrows!(plot, points, direcs,
    linewidth  = axesbody,
    arrowsize  = axeshead,
    arrowcolor = axescolor,
    linecolor  = axescolor,
  )

  # plot axes names
  position = Tuple.(direcs)
  Makie.text!(plot, axesnames,
    position = position,
  )

  # plot samples
  points = [Makie.Point{d}(v) for v in eachrow(F)]
  Makie.scatter!(plot, points,
    markersize = dotsize,
    color = dotcolor,
  )

  if showdotnames
    # plot sample names
    position = Tuple.(points)
    Makie.text!(plot, string.(1:n),
      position = position,
    )
  end
end

end
