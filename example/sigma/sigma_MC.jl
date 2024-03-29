using ElectronLiquid
using CompositeGrids

include("../input.jl")

# mission = :Z
# mission = :K
mission = ARGS[1]
println("mission: ", mission)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock, spin=spin)
    println(UEG.short(para))
    # para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock, spin=1)
    kF = para.kF

    if mission == "Z"
        ######### calcualte Z factor ######################
        kgrid = [kF,]
        # ngrid = [-1, 0, 1]
        ngrid = [0]
    elseif mission == "K"
        ######### calculate K dependence #####################
        Nk, korder = 4, 4
        minK = 0.2kF
        # kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.2kF], [kF,], Nk, minK, korder).grid
        # kgrid = kF .+ [-0.1, -0.05, -0.03, -0.01, -0.005, -0.001, 0, 0.001, 0.005, 0.01, 0.03, 0.05, 0.1] * kF
        kgrid = kF .+ [-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1] * kF
        ngrid = [0,]
    else
        error("unknown mission")
    end

    if isLayered2D
        filename = "./data$(dim)d/data_$(mission)_layered2d.jld2"
    else
        # filename = "./data$(dim)d/data$(dim)d_$(mission)_test.jld2"
        filename = "./data$(dim)d/data_Z_test.jld2"
    end
    # filename = "data$(dim)_$(mission).jld2"

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

    sigma, result = Sigma.MC(para; kgrid=kgrid, ngrid=ngrid, spinPolarPara=spinPolarPara,
        neval=neval, filename=filename, partition=partition, reweight_goal=reweight_goal,
        # isLayered2D=isLayered2D,
        diagtype=diagGenerate)
    # isClib=false)
end
