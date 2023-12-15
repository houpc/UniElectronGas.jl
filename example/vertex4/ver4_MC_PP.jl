using ElectronLiquid, FeynmanDiagram, JLD2

dim = 3
# rs = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
# rs = [1.0, 2.0, 3.0, 4.0, 5.0]
rs = [0.5,]
# Fs = -[0.223, 0.380, 0.516, 0.639, 0.752]
# Fs = -[0.223,]
# mass2 = [1e-3,]
# mass2 = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
# mass2 = [4.0, 5.0, 6.0]
mass2 = [4.0,]
# Fs = [-0.0,]
Fs = -0.0 .* rs
beta = [100.0]
order = [3,]
# ell = [0, 1]
ell = 0
neval = 1e6
# neval = 1e8
# isDynamic = true
isDynamic = false
isFock = false

# for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
for (irs, _mass2, _beta, _order) in Iterators.product([i for i in 1:length(rs)], mass2, beta, order)
    _F = Fs[irs]
    _rs = rs[irs]
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    println(UEG.short(para))
    # filename = "data_ver4PP.jld2"
    # filename = "data_ver4PP_profile.jld2"
    # filename = "data_ver4PP_new.jld2"
    filename = "data_ver4PP_parqAD.jld2"

    partition = UEG.partition(_order)
    # partition = partition[1:end] # exclude order 1
    # partition = [(1, 0, 2), (1, 1, 1), (1, 2, 0), (2, 0, 1), (2, 1, 0), (3, 0, 0)]
    # partition = [(1, 0, 0), (2, 0, 0), (2, 1, 0), (3, 0, 0)]
    println(partition)
    neighbor = UEG.neighbor(partition)
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in partition
        order == 1 && sOrder > 0 && continue
        push!(reweight_goal, 2.0^(2order + sOrder + vOrder - 2))
    end
    push!(reweight_goal, 4.0)
    # reweight_goal = [1.0, 4.0, 8.0, 16.0, 4.0]
    println(reweight_goal)

    ver4, result = Ver4.MC_PP_ParquetAD(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition, filter=[NoHartree,])
    # ver4, result = Ver4.MC_PP(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition, filter=[NoHartree, NoBubble])
    # ver4, result = Ver4.MC_PP(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition, filter=[NoHartree,])
    # ver4, result = Ver4.MC_PH(para; l=ell, neval=neval, filename=filename, partition=partition, filter=[NoHartree])
    # jldopen("test.jld2", "w") do f
    #     key = "test"
    #     if haskey(f, key)
    #         @warn("replacing existing data for $key")
    #         delete!(f, key)
    #     end
    #     fake = Dict((1, 0, 0) => [1, 2, 3], (2, 0, 0) => [4, 5, 6], (3, 0, 0) => "3")
    #     f["test"] = fake
    #     println(ver4)
    #     for k in keys(ver4)
    #         println(k)
    #         println(ver4[k])
    #         f["$k"] = ver4[k]
    #     end
    # end

end
