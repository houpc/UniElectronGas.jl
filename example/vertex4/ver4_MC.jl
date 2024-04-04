using ElectronLiquid, FeynmanDiagram
import FeynmanDiagram.FrontEnds: PHr, PHEr, PPr, Alli

# dim = 2
# rs = [0.5,]
# mass2 = [1.0,]
# Fs = [-0.0,]
# beta = [25.0]
# order = [4,]
ell = [0,]
channels = [PHr, PHEr, PPr, Alli]
# channels = [Alli,]
# channels = [PHr,]
# channels = [PHEr,]
# channels = [PPr,]
# neval = 2e7
# neval = 1e6
# isDynamic = false
# isFock = false

include("../input.jl")

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    println(UEG.short(para))
    filename = "data$(dim)d_ver4PH_test.jld2"

    # _partition = UEG.partition(_order, offset=0)
    _partition = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 0, 0)]
    # _partition = [(3, 0, 0), (4, 0, 0)]
    partition = Vector{Tuple{Int64,Int64,Int64}}()
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in _partition
        order == 0 && sOrder > 0 && continue
        push!(partition, (order, sOrder, vOrder))
        push!(reweight_goal, 2.0^(2order + sOrder + vOrder))
    end
    push!(reweight_goal, 1.0)
    neighbor = UEG.neighbor(partition)

    # ver4, result = Ver4.MC_PH_AD(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition,
    # generate_type=diagGenerate, channels=channels)
    ver4, result = Ver4.MC_PH_Clib(para; l=ell, neval=neval, filename=filename, reweight_goal=reweight_goal, partition=partition,
        generate_type=diagGenerate)
end
