using ElectronLiquid, FeynmanDiagram

dim = 3
# rs = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
# rs = [1.0, 2.0, 3.0, 4.0, 5.0]
rs = [2.0,]
# Fs = -[0.223, 0.380, 0.516, 0.639, 0.752]
# Fs = -[0.223,]
# mass2 = [1e-3,]
# mass2 = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
# mass2 = [1.0, 2.0, 3.0, 4.0, 5.0]
mass2 = [2.5, 3.5]
# mass2 = [3.6, 3.8]
Fs = [-0.0,]
# Fs = -0.2 .* rs
beta = [100.0]
order = [4,]
# ell = [0, 1]
ell = 0
neval = 2e7
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
    filename = "data_ver4PP_profile.jld2"

    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in partition
        order == 1 && sOrder > 0 && continue
        push!(reweight_goal, 2.0^(2order + sOrder + vOrder - 2))
    end
    push!(reweight_goal, 4.0)

    # ver4, result = Ver4.MC_PP(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition, filter=[NoHartree, NoBubble])
    ver4, result = Ver4.MC_PP(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition, filter=[NoHartree,])
    # ver4, result = Ver4.MC_PH(para; l=ell, neval=neval, filename=filename, partition=partition, filter=[NoHartree])
end
