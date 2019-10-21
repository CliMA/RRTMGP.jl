ENV["DATADEPS_ALWAYS_ACCEPT"] = true

using Test
using DataDeps
using JRRTMGP

@testset "DataDep" begin
  testfolder = JRRTMGP.data_folder_rrtmgp_test_val()
  @test testfolder isa String
end

