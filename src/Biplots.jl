module Biplots

using LinearAlgebra
using Statistics

import Makie

"""
    biplot(X)

Biplot of design matrix `X`. See https://en.wikipedia.org/wiki/Biplot.

* `α` - Shape parameter `α ∈ [0,1]`

## References

* Gabriel, K. 1971. [The Biplot Graphic Display of Matrices
  with Application to Principal Component Analysis]
  (https://www.jstor.org/stable/2334381)

* Aitchison, J. & Greenacre, M. 2002. [Biplots of Compositional data]
  (https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9876.00275)
"""
@Makie.recipe(Biplot, X) do scene
  Makie.Attributes(;
    # generic attributes
    markersize = Makie.theme(scene, :markersize),

    # custom attributes
    d = 2,   # number of dimensions
    α = 1.0, # shape parameter
  )
end

function Makie.plot!(plot::Biplot{<:Tuple{AbstractMatrix}})
  # retrieve parameters
  X = plot[:X][]
  d = plot[:d][]
  α = plot[:α][]

  # sanity checks
  @assert d ∈ [2,3] "d must be 2 or 3"
  @assert 0 ≤ α ≤ 1 "α must be in [0,1]"

  # number of samples
  n = size(X, 1)

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
  Makie.arrows!(plot, points, direcs)

  # plot samples
  points = [Makie.Point{d}(v) for v in eachrow(F)]
  Makie.scatter!(plot, points)
end

end
