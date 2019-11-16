using Test
using DataDeps
using RRTMGP

@testset "DataDep" begin
  testfolder = RRTMGP.data_folder_rrtmgp_test_val()
  @test testfolder isa String
end

