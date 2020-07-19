cd(@__DIR__)

using Pkg; Pkg.activate(".."); Pkg.instantiate()

using CSV, DataFrames, DelimitedFiles, ProgressMeter, Pipe

#---

dsts = readdlm("../data/levenshteinLanguageDistances.csv", ',')

languages = string.(sort(unique(vcat(dsts[1,2:end], dsts[2:end,1]))))

nonProto = filter(x -> !occursin("PROTO", x), languages)

nPairs = length(nonProto)*(length(nonProto)-1) ÷ 2


dstArray = zeros(length(languages), length(languages))
for j in 1:length(languages)
    for i in 2:length(languages)
        if j>=i
            dstArray[i-1,j] = dsts[i,j]
        end
    end
end

prInd = Tuple.(findall(dstArray .<= 0.7))
filter!(x -> x[1] < x[2], prInd)

probablyRelated = @pipe [
          [languages[x[1]], languages[x[2]]]
          for
          x in prInd if languages[x[1]] ∈ nonProto && languages[x[2]] ∈ nonProto
      ] |>
      hcat(_...) |>
      permutedims


writedlm("../data/probablyRelated.csv", probablyRelated, ',')

#---
d = CSV.File("../data/asjp19wide.csv") |> DataFrame!

prI = indexin(probablyRelated, d.longname)

wordPairs = Tuple{String, String}[]

@showprogress for i in 1:size(prI,1)
    j,k = prI[i,:]
    wl1 = Vector(d[j,2:end])
    wl2 = Vector(d[k,2:end])
    for (w1,w2) in zip(wl1, wl2)
        if !ismissing(w1) && !ismissing(w2)
            for x in split(w1,"-")
                for y in split(w2,"-")
                    push!(wordPairs, (x,y))
                end
            end
        end
    end
end

@time wordPairs = permutedims(hcat(collect.(wordPairs)...))

writedlm("../data/potentialCognates.csv", wordPairs, ',')
