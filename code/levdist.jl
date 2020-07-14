cd(@__DIR__)

using Pkg; Pkg.activate(".."); Pkg.instantiate()

#using Revise

include("SequenceAlignment/SequenceAlignment.jl")

using CSV, DataFrames, LinearAlgebra, StatsBase
import .SequenceAlignment:ldn

# extend ldn to missing
function ldn(w1::Union{Missing,String}, w2::Union{Missing,String})
    if ismissing(w1) || ismissing(w2)
        return missing
    else
        return SequenceAlignment.ldn(w1,w2)
    end
end

function ldnStar(w1::Union{Missing,String},w2::Union{Missing,String})
    if ismissing(w1) || ismissing(w2)
        return missing
    end
    minimum([ldn(x,y) for x in split(w1,"-") for y in split(w2,"-")])
end



const minSim = -sqrt(40)
const maxSim = (log(40*(39)+1)-1)*sqrt(40)


const d = CSV.File("../data/asjp19wide.csv") |> DataFrame!

const languages = string.(d.longname)

const lInd = indexmap(languages)

const concepts = names(d)[2:end]

const asjpMtx = Matrix(d[:,2:end])


function dercLDN(
    l1::String,
    l2::String
    # asjpMtx::Matrix{Union{Missing, String}} = asjpMtx,
    # maxSim::Float64 = maxSim,
    # minSim::Float64 = minSim,
)
    i, j = lInd[l1], lInd[l2]
    defined1 = findall(.!ismissing.(asjpMtx[i, :]))
    defined2 = findall(.!ismissing.(asjpMtx[j, :]))
    definedBoth = intersect(defined1, defined2)
    nBoth = length(definedBoth)
    dg = Vector{Float64}(undef, nBoth)
    for (k, c) in enumerate(definedBoth)
        dg[k] = ldnStar(asjpMtx[i, c], asjpMtx[j, c])
    end
    offDg = Vector{Float64}(undef, length(defined1) * length(defined2) - nBoth)
    counter = 1
    for k1 in defined1
        w1 = asjpMtx[i, k1]
        for k2 in defined2
            if k1 != k2
                w2 = asjpMtx[j, k2]
                offDg[counter] = ldnStar(w1, w2)
                counter += 1
            end
        end
    end
    ranks = Vector{Float64}(undef, nBoth)
    for k = 1:nBoth
        x = dg[k]
        ranks[k] = geomean(1 .+ (sum(offDg .< x):sum(offDg .<= x)))
    end
    stc = mean(-log.(ranks ./ (1 + length(offDg))))
    sim = (stc - 1) * sqrt(length(dg))
    (maxSim - sim) / (maxSim - minSim)
end


ger = "IE.GERMANIC.STANDARD_GERMAN"

eng = "IE.GERMANIC.ENGLISH"

@show dercLDN(ger, eng)

@time gerDsts = [dercLDN(ger, l) for l in languages]