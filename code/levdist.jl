cd(@__DIR__)


using Pkg; Pkg.activate(".."); Pkg.instantiate()

const nthreads = 50
using Distributed
using SharedArrays

addprocs(nthreads)

@everywhere mean(a) = sum(a)/length(a)

@everywhere geomean(a) = exp(mean(log.(a)))


@everywhere include("SequenceAlignment/SequenceAlignment.jl")

using CSV, DataFrames, LinearAlgebra
using SharedArrays
@everywhere using .SequenceAlignment
#@everywhere using StatsBase

@everywhere function ldnStar(w1::Union{Missing,String},w2::Union{Missing,String})
    if ismissing(w1) || ismissing(w2)
        return missing
    end
    v1 = split(w1,"-")
    v2 = split(w2,"-")
    scores = Vector{Float64}(undef, length(v1)*length(v2))
    counter = 1
    for x in v1
        for y in v2
            @inbounds scores[counter] = ldn(x,y)
            counter += 1
        end
    end
    minimum(scores)
end





@everywhere begin
    using CSV
    using DataFrames
    d = CSV.File("../data/asjp19wide.csv") |> DataFrame
    minSim = -sqrt(40)
    maxSim = (log(40*(39)+1)-1)*sqrt(40)
end;

languages = string.(d.longname);

# concepts = names(d)[2:end]


@everywhere asjpMtx = Matrix(d[:,2:end]);



@everywhere function dercLDN(
    i::Int,
    j::Int,
    asjpMtx::Matrix{Union{Missing, String}} = asjpMtx,
    maxSim::Float64 = maxSim,
    minSim::Float64 = minSim,
)
    defined1 = findall(.!ismissing.(asjpMtx[i, :]))
    defined2 = findall(.!ismissing.(asjpMtx[j, :]))
    definedBoth = intersect(defined1, defined2)
    nBoth = length(definedBoth)
    dg = Vector{Float64}(undef, nBoth)
    for (k, c) in enumerate(definedBoth)
        dg[k] = ldnStar(asjpMtx[i, c], asjpMtx[j, c])
    end
    nOffD = length(defined1) * length(defined2) - nBoth
    offDg = Vector{Float64}(undef, nOffD)
    counter = 1
    for k1 in defined1
        w1 = asjpMtx[i, k1]
        for k2 in defined2
            if k1 != k2
                w2 = asjpMtx[j, k2]
                @inbounds offDg[counter] = ldnStar(w1, w2)
                counter += 1
            end
        end
    end
    ranks = Vector{Float64}(undef, nBoth)
    for k = 1:nBoth
        @inbounds x = dg[k]
        @inbounds ranks[k] = geomean(1 .+ (sum(offDg .< x):sum(offDg .<= x)))
    end
    stc = mean(-log.(ranks ./ (1 + nOffD)))
    sim = (stc - 1) * sqrt(nBoth)
    (maxSim - sim) / (maxSim - minSim)
end


# ger = "IE.GERMANIC.STANDARD_GERMAN"
#
# eng = "IE.GERMANIC.ENGLISH"
#
# @show dercLDN(ger, eng)
#
# @time gerDsts = [dercLDN(ger, l) for l in languages]
#
# @time begin
#     gerDistsT = Vector{Union{Missing, Float64}}(undef, length(languages))
#     Threads.@threads for i in 1:length(gerDistsT)
#         gerDistsT[i] = dercLDN(ger, languages[i])
#     end
# end


lpairs = [
    [i,j] for (i, l1) in enumerate(languages)
    for (j, l2) in enumerate(languages) if i < j
];



@everywhere function f(p)
    i,j = p
    dercLDN(i,j)
end

#---

ldists = pmap(f, lpairs);


levDistsLong = DataFrame(l1=languages[hcat(lpairs...)[1,:]],
                         l2=languages[hcat(lpairs...)[2,:]], ldist=ldists)

levDists = unstack(levDistsLong, :l1, :l2, :ldist)

CSV.write("../data/levenshteinLanguageDistances.csv", levDists)
