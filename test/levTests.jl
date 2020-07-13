cd(@__DIR__)

using CSV, Test, DataFrames

include("../code/alignment.jl")

using .SequenceAlignment

d = CSV.File("levenshteinTests.csv") |> DataFrame!

@testset "Levenshtein Distance" begin
    for rw in eachrow(d)
        w1, w2, l = rw
        @test levenshteinDistance(w1, w2) == l
    end
end
