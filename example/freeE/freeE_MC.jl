using ElectronLiquid
using JLD2

include("../input.jl")

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    println(UEG.short(para))

    if isLayered2D
        filename = "data_freeE_layered2d.jld2"
    else
        filename = "data_freeE.jld2"
    end

    # freeE, result = FreeEnergy.MC(para; neval=neval, filename=filename, spinPolarPara=spinPolarPara)
    freeE, result = FreeEnergy.MC(para;
        isLayered2D=isLayered2D,
        neval=neval, filename=filename,
        spinPolarPara=spinPolarPara)

    # freeE, result = FreeEnergy.MC(para; neval=neval, filename=filename, spinPolarPara=spinPolarPara)
    # freeE, result = FreeEnergy.MC(para; neval=neval, filename=filename,
    #     spinPolarPara=spinPolarPara, partition=[(0, 0, 0)])
end