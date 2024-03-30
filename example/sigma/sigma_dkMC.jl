using ElectronLiquid
using CompositeGrids

include("../input.jl")

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock, spin=spin)
    println(UEG.short(para))
    # para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock, spin=1)
    kgrid = [para.kF]
    ngrid = [0,]

    partition = UEG.partition(_order)
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in partition
        reweight_factor = 2.0^(2order + 2sOrder + vOrder - 2)
        if (order, sOrder, vOrder) == (1, 0, 0)
            reweight_factor = 6.0
        end
        push!(reweight_goal, reweight_factor)
    end
    push!(reweight_goal, 4.0)

    filename = mission == "Z" ? sigma_z_filename : sigma_k_filename

    # sigma, result = Sigma.MC_dk(para; kgrid=kgrid, ngrid=ngrid,
    sigma, result = Sigma.MC_dk_Clib(para; kgrid=kgrid, ngrid=ngrid,
        neval=neval, filename=filename, partition=partition, reweight_goal=reweight_goal,
        diagtype=diagGenerate, isLayered2D=isLayered2D)
end
