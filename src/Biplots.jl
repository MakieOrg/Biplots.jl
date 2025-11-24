# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Biplots

using Tables
using LinearAlgebra
using Statistics
using Printf

import Makie

# -----------------------
# TYPES OF NORMALIZATION
# -----------------------

center(X) = X .- mean(X, dims=1)

function logcenter(X)
  L  = log.(X)
  μn = mean(L, dims=1)
  μp = mean(L, dims=2)
  μ  = mean(L)
  L .- μn .- μp .+ μ
end

formtransform(X)  = (Z=center(X),    α=1.0, κ=√(size(X,1)-1))
covtransform(X)   = (Z=center(X),    α=0.0, κ=1.0           )
rformtransform(X) = (Z=logcenter(X), α=1.0, κ=1.0           )
rcovtransform(X)  = (Z=logcenter(X), α=0.0, κ=√size(X,2)    )

biplotof = Dict(
  :form => formtransform,
  :cov  => covtransform,
  :rform => rformtransform,
  :rcov => rcovtransform
)

"""
    factors(table; kind=:form, dim=2)

Return `F`, `G/k` and `σ²` where `F` and `G` are the
factors (i.e. matrices with rank `dim`) of the design
matrix `Z` obtained with a given `kind` of transformation.

The scalar `κ` scales the column vectors for plotting
purposes and the variance explained `σ²` is also returned.
"""
function factors(table; kind=:form, dim=2)
  # sanity checks
  @assert dim ∈ [2,3] "dim must be 2 or 3"
  @assert kind ∈ [:form,:cov,:rform,:rcov] "$kind is not a valid kind of biplot"

  # original design matrix
  X = Tables.matrix(table)

  # matrix transformation
  Z, α, κ = biplotof[kind](X)

  # singular value decomposition
  U, σ, V = svd(Z)

  # matrix factors Z ≈ F*G'
  F = U[:,1:dim] .* (σ[1:dim] .^ α)'
  G = V[:,1:dim] .* (σ[1:dim] .^ (1 - α))'

  # variance explained
  σ² = σ[1:dim] .^ 2 / sum(σ .^ 2)

  # return scaled factors
  F, G/κ, σ²
end

"""
    biplot(table; kind=:form, [aesthetics...])

Biplot of `table` of given `kind` in 2-dimensional space.

There are four kinds of biplots:

* `:form`  - standard biplot with shape parameter `α = 1`
* `:cov`   - standard biplot with shape parameter `α = 0`
* `:rform` - relative variation biplot with `α = 1`
* `:rcov`  - relative variation biplot with `α = 0`

# Biplot attributes

* `kind` - kind of biplot (`:form`, `:cov`, `:rform` or `:rcov`)

# Aesthetics attributes

* `axescolormap` - colormap of principal axes (default to theme colormap)
* `dotcolormap`  - colormap of sample dots (default to theme colormap)
* `fontsize`     - font size in axes and dots (default to `12`)
* `axesbody`     - size of principal axes' body (default to `2`)
* `axeshead`     - size of principal axes' head (default to `8`)
* `axescolor`    - color of principal axes (default to `:black`)
* `axeslabel`    - names of principal axes (default to `x1,x2,...`)
* `dotsize`      - size of sample dots (default to `4`)
* `dotcolor`     - color of sample dots (default to `:black`)
* `dotlabel`     - names of sample dots (default to `1:nobs`)
* `showdots`     - show names of dots (default to `true`)
* `showlinks`    - show links between principal axes (default to`false`)

See https://en.wikipedia.org/wiki/Biplot.

## References

* Gabriel, K. 1971. [The Biplot Graphic Display of Matrices
  with Application to Principal Component Analysis]
  (https://www.jstor.org/stable/2334381)

* Aitchison, J. & Greenacre, M. 2002. [Biplots of Compositional data]
  (https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/1467-9876.00275)
"""
@Makie.recipe(Biplot, table) do scene
  Makie.Attributes(;
    # biplot attributes
    kind = :form,

    # aesthetic attributes
    axescolormap = Makie.theme(scene, :colormap),
    dotcolormap  = Makie.theme(scene, :colormap),
    fontsize     = 12,
    axesbody     = 2,
    axeshead     = 8,
    axescolor    = :black,
    dotsize      = 4,
    dotcolor     = :black,
    dotlabel     = nothing,
    showdots     = true,
    showlinks    = false,
  )
end

function Makie.plot!(plot::Biplot{<:Tuple{Any}})
  # biplot attributes
  table = plot[:table][]
  kind  = plot[:kind][]

  # aesthetics attributes
  axescolormap = plot[:axescolormap][]
  dotcolormap  = plot[:dotcolormap][]
  fontsize     = plot[:fontsize][]
  axesbody     = plot[:axesbody][]
  axeshead     = plot[:axeshead][]
  axescolor    = plot[:axescolor][]
  dotsize      = plot[:dotsize][]
  dotcolor     = plot[:dotcolor][]
  dotlabel     = plot[:dotlabel][]
  showdots     = plot[:showdots][]
  showlinks    = plot[:showlinks][]

  # dimension of the biplot
  dim = 2

  # biplot factors
  F, G, σ² = factors(table, kind=kind, dim=dim)

  # n samples x p variables
  n = size(F, 1)
  p = size(G, 1)

  # defaults aesthetics
  if axescolor isa AbstractVector{<:Number}
    min, max = extrema(axescolor)
    axesscale(x) = x / (max-min) - min / (max - min)
    axescolor = Makie.cgrad(axescolormap)[axesscale.(axescolor)]
  end
  if isnothing(dotlabel)
    dotlabel = string.(1:n)
  end
  if dotcolor isa AbstractVector{<:Number}
    min, max = extrema(dotcolor)
    dotscale(x) = x / (max-min) - min / (max - min)
    dotcolor = Makie.cgrad(dotcolormap)[dotscale.(dotcolor)]
  end

  # plot principal axes
  origin = fill(Makie.Point(ntuple(i->0.0, dim)), p)
  direcs = [Makie.Vec{dim}(v) for v in eachrow(G)]
  Makie.arrows2d!(plot, origin, direcs,
    shaftwidth = axesbody,
    tipwidth = axeshead,
    color = axescolor
  )

  # plot axes names
  cols = Tables.columns(table)
  names = [string(n) for n in Tables.columnnames(cols)]
  Makie.annotation!(plot, direcs,
    text = names,
    textcolor = axescolor,
    fontsize = fontsize
  )

  # plot links between axes
  if showlinks
    links = Makie.Vec{dim}[]
    for i in 1:p
      for j in i+1:p
        push!(links, direcs[i])
        push!(links, direcs[j])
      end
    end
    position = Tuple.(links)
    Makie.linesegments!(plot, position,
      color = :gray,
      linestyle = :dash,
      linewidth = 0.5,
    )
  end

  # plot samples
  points = [Makie.Point{dim}(v) for v in eachrow(F)]
  Makie.scatter!(plot, points,
    markersize = dotsize,
    color = dotcolor,
  )

  # plot samples names
  if showdots
    Makie.annotation!(plot, points,
      text = dotlabel,
      textcolor = dotcolor,
      fontsize = fontsize
    )
  end

  # plot variance explained
  minpos = fill(Inf, dim)
  for i in 1:dim
    for direc in direcs
      if direc[i] < minpos[i]
        minpos[i] = direc[i]
      end
    end
    for point in points
      if point[i] < minpos[i]
        minpos[i] = point[i]
      end
    end
  end
  textdim = [(@sprintf "Dim %d (%.01f " i 100*v)*"%)" for (i, v) in enumerate(σ²)]
  Makie.text!(plot, join(textdim, "\n"),
    position = Tuple(minpos),
  )
end

end
