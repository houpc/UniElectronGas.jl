# using ElectronLiquid, FeynmanDiagram, 
using JLD2, Measurements

jldopen("test.jld2", "w") do f
    key = "test"
    if haskey(f, key)
        @warn("replacing existing data for $key")
        delete!(f, key)
    end
    f["test"] = [measurement(8.7e10, 1.5e11)+measurement(1.3e11, 1.0e11)*1im; measurement(8.7e10, 1.5e11)+measurement(1.3e11, 1.0e11)*1im;;;]
end

jldopen("test.jld2", "r") do f
    println(keys(f))
    for k in keys(f)
        println(k)
        println(f[k])
    end
end
