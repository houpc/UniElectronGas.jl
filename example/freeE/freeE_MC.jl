using ElectronLiquid
using JLD2

include("../input.jl")

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, spin=spin, isFock=isFock)
    println(UEG.short(para))

    _partition = UEG.partition(_order, offset=0)
    partition = Vector{Tuple{Int,Int,Int}}()
    for p in _partition
        isLayered2D && p[3] > 0 && continue
        push!(partition, p)
    end
    println("partition: ", partition)

    if isLayered2D
        filename = "data_freeE_layered2d_hpc.jld2"
    else
        filename = "data$(dim)d_freeE.jld2"
    end

    # freeE, result = FreeEnergy.MC(para;
    #     isLayered2D=isLayered2D, partition=partition,
    #     neval=neval, filename=filename)

    freeE, result = FreeEnergy.MC_Clib(para;
        isLayered2D=isLayered2D, partition=partition,
        neval=neval, filename=filename)
end