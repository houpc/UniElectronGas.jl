using ElectronLiquid
using JLD2

dim = 2
rs = [1.0]
mass2 = [1.0,]
# mass2 = [0.5, 1.0, 2.0, 4.0]
# mass2 = [0.001]
Fs = [-0.0,]
# beta = [25.0, 40.0, 80.0]
beta = [25.0]
order = [0,]
neval = 1e6
isDynamic = false # dynamically screened Coulomb interaction or not
isFock = false # Fock renormalized Green's function or not
spinPolarPara = 0.0 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
isLayered2D = true

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