using Biplots
using CairoMakie
using ReferenceTests
using Test

@testset "Biplots.jl" begin
  # data matrix (22 paintings x 6 colors)
  data = [
    0.125 0.243 0.153 0.031 0.181 0.266
    0.143 0.224 0.111 0.051 0.159 0.313
    0.147 0.231 0.058 0.129 0.133 0.303
    0.164 0.209 0.120 0.047 0.178 0.282
    0.197 0.151 0.132 0.033 0.188 0.299
    0.157 0.256 0.072 0.116 0.153 0.246
    0.153 0.232 0.101 0.062 0.170 0.282
    0.115 0.249 0.176 0.025 0.176 0.259
    0.178 0.167 0.048 0.143 0.118 0.347
    0.164 0.183 0.158 0.027 0.186 0.281
    0.175 0.211 0.070 0.104 0.157 0.283
    0.168 0.192 0.120 0.044 0.171 0.305
    0.155 0.251 0.091 0.085 0.161 0.257
    0.126 0.273 0.045 0.156 0.131 0.269
    0.199 0.170 0.080 0.076 0.158 0.318
    0.163 0.196 0.107 0.054 0.144 0.335
    0.136 0.185 0.162 0.020 0.193 0.304
    0.184 0.152 0.110 0.039 0.165 0.350
    0.169 0.207 0.111 0.057 0.156 0.300
    0.146 0.240 0.141 0.038 0.184 0.250
    0.200 0.172 0.059 0.120 0.136 0.313
    0.135 0.225 0.217 0.019 0.187 0.217
  ]

  # variable names
  names = [:Black,:White,:Blue,:Red,:Yellow,:Other]

  # choose any Tables.jl table
  table = (; zip(names, eachcol(data))...)

  fig, ax = biplot(table, kind = :form, dotcolor = table.Red)
  ax.aspect = DataAspect()
  @test_reference "data/form.png" fig

  fig, ax = biplot(table, kind = :cov, dotcolor = table.Red)
  ax.aspect = DataAspect()
  @test_reference "data/cov.png" fig

  fig, ax = biplot(table, kind = :rform, dotcolor = table.Red)
  ax.aspect = DataAspect()
  @test_reference "data/rform.png" fig

  fig, ax = biplot(table, kind = :rcov, dotcolor = table.Red)
  ax.aspect = DataAspect()
  @test_reference "data/rcov.png" fig
end
