using ElectronLiquid, FeynmanDiagram

dim = 3
rs = [1.0,]
mass2 = [1e-3,]
Fs = [-0.0,]
beta = [25.0]
order = [2,]
# ell = [0, 1]
ell = 0
# neval = 2e7
neval = 1e6
isDynamic = true
isFock = false

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    println(UEG.short(para))
    filename = "data_ver4PP.jld2"

    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in partition
        order == 1 && sOrder > 0 && continue
        push!(reweight_goal, 2.0^(2order + sOrder + vOrder - 2))
    end
    push!(reweight_goal, 4.0)

    ver4, result = Ver4.MC_PP(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition, filter=[NoHartree, NoBubble])
    # ver4, result = Ver4.MC_PH(para; l=ell, neval=neval, filename=filename, partition=partition, filter=[NoHartree])
end
