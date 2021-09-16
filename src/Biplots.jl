module Biplots

using LinearAlgebra
using Statistics

import Makie

"""
    biplot(X)

Biplot of design matrix `X` with various options:

# Biplot attributes

* `d` - number of dimensions `d ∈ {2,3}`
* `α` - shape parameter `α ∈ [0,1]`

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
    d = 2,
    α = 1.0,

    # aesthetic attributes
    axesbody  = nothing,
    axeshead  = nothing,
    axescolor = :black,
    axesnames = nothing,
    dotsize   = nothing,
    dotcolor  = :black,
  )
end

function Makie.plot!(plot::Biplot{<:Tuple{AbstractMatrix}})
  # retrieve parameters
  X = plot[:X][]
  d = plot[:d][]
  α = plot[:α][]

  # retrieve options
  axesbody  = plot[:axesbody][]
  axeshead  = plot[:axeshead][]
  axescolor = plot[:axescolor][]
  axesnames = plot[:axesnames][]
  dotsize   = plot[:dotsize][]
  dotcolor  = plot[:dotcolor][]

  # size of design matrix
  n, m = size(X)

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

  # sanity checks
  @assert d ∈ [2,3] "d must be 2 or 3"
  @assert 0 ≤ α ≤ 1 "α must be in [0,1]"
  @assert length(axesnames) == m "axesnames must have length $m"

  # center matrix
  Z = X .- mean(X, dims=1)

  # singular value decomposition
  U, σ, V = svd(Z)

  # variance explained
  v = σ[1:d] .^ 2 / sum(σ .^ 2)

  # matrix factors X ≈ F*G'
  F = U[:,1:d] .* (σ[1:d] .^ α)'
  G = V[:,1:d] .* (σ[1:d] .^ (1 - α))'

  # plot principal axes
  points = fill(Makie.Point(ntuple(i->0., d)), n)
  direcs = [Makie.Vec{d}(v) for v in eachrow(G)]
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
end

end
