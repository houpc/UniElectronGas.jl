using ElectronLiquid, FeynmanDiagram

# dim = 2
# rs = [0.5,]
# mass2 = [1.0,]
# Fs = [-0.0,]
# beta = [25.0]
# order = [4,]
ell = [0, 1]
# neval = 2e7
# neval = 1e6
# isDynamic = false
# isFock = false

include("../input.jl")

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    println(UEG.short(para))
    filename = "data_ver4PH.jld2"

    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in partition
        order == 1 && sOrder > 0 && continue
        push!(reweight_goal, 2.0^(2order + sOrder + vOrder - 2))
    end
    push!(reweight_goal, 4.0)

    ver4, result = Ver4.MC_PH(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition, filter=[NoHartree])
    # ver4, result = Ver4.MC_PH(para; l=ell, neval=neval, filename=filename, partition=partition, filter=[NoHartree])
end
