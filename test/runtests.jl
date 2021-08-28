using GridapGeosciences
using Test

@testset "GridapGeosciences" begin
  @time @testset "CubedSphereDiscreteModelsTests" begin include("CubedSphereDiscreteModelsTests.jl") end
  @time @testset "DarcyCubedSphereTests" begin include("DarcyCubedSphereTests.jl") end
  @time @testset "WeakDivPerpTests" begin include("WeakDivPerpTests.jl") end
  @time @testset "LaplaceBeltramiCubedSphereTests" begin include("LaplaceBeltramiCubedSphereTests.jl") end
  @time @testset "WaveEquationCubedSphereTests" begin include("WaveEquationCubedSphereTests.jl") end
  @time @testset "ThetaMethodFullNewtonNSWECubedSphereTests" begin include("ThetaMethodFullNewtonNSWECubedSphereTests.jl") end
end
