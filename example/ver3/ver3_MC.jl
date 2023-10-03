
include("ver3_static.jl")

dim = 3

rs = [1.0,]
# mass2 = [1.0, 2.0, 3.0, 4.0, 5.0]
mass2 = [5.0,]
# mass2 = [1.0,]
# mass2 = [6.0, 8.0, 10.0, 12.0, 14.0]
# mass2 = [10.5, 11.0]
# mass2 = [1e-3,]
Fs = [-0.0,]
# Fs = -[0.223, 0.380, 0.516, 0.639, 0.752]
# Fs = -[0.223,]
beta = [100.0]
order = [4,]
ell = 0
neval = 1e7
# isDynamic = true
isDynamic = false
isFock = false

const filename = "data_ver3.jld2"

# anglegrid = [0.0, 0.25π, 0.5π, 0.75π, π]
anglegrid = [π,]

for (irs, _mass2, _beta, _order) in Iterators.product([i for i in 1:length(rs)], mass2, beta, order)
    _F = Fs[irs]
    _rs = rs[irs]
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    println(UEG.short(para))

    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in partition
        # order == 1 && sOrder > 0 && continue
        push!(reweight_goal, 2.0^(2order + sOrder + vOrder - 2))
    end
    push!(reweight_goal, 4.0)
    channel = [PHr, PHEr, PPr]
    diagram = Ver3.diagram(para, partition; filter=[NoHartree, Proper])
    # diagram = Ver3.diagram(para, partition; filter=[NoHartree, NoBubble, Proper])

    ver3, result = ver3_static(para, diagram, anglegrid; neval=neval, reweight_goal=reweight_goal, filename=filename)
    # ver3, result = ver3_static(para, diagram, anglegrid; neval=neval, filename=filename)

end
