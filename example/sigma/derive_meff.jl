using UniElectronGas, ElectronLiquid
using JLD2, DelimitedFiles

include("../input.jl")

# # dim = 3
# dim = 2
# spin = 1
# # spin = 2
# # rs = [0.5, 1.0, 2.0]
# rs = [1.0]
# # rs = [0.1, 0.3]
# # mass2 = [0.5, 1.0, 1.5, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
# # mass2 = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
# mass2 = [1.75]
# Fs = [-0.0,]
# # beta = [25.0,]
# beta = [40.0,]
# order = [4,]
# # order = [5,]
# isDynamic = false
### dSigma/dk = (Sigma[kF_label+idx_dk] - Sigma[kF_label-idx_dk]) / (kgrid[kF_label+idx_dk] - kgrid[kF_label-idx_dk])

# const parafilename = "para_wn_1minus0.csv"
const filename = "./data$(dim)d/data$(dim)d_K.jld2"
# const filename = "./data$(dim)d/data$(dim)d_K_o5.jld2"
const savefilename = spin == 2 ? "meff_$(dim)d.dat" : "meff_$(dim)d_spin$spin.dat"

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")
    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, spin=spin)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                println(UEG.paraid(para))

                ngrid, kgrid, rSw_k, iSw_k = UniElectronGas.getSigma(para, filename)
                meff, fit_p = UniElectronGas.getMeff(para, rSw_k, kgrid)
                # meff = UniElectronGas.getMeff(para, rSw_k, kgrid)
                println("m* / m = ", meff)
                push!(results, Any[_rs, _beta, _mass2, _order, meff...])
                writedlm("fit.dat", fit_p)
            end
        end
    end
    if isSave
        open(savefilename, "a+") do io
            writedlm(io, results)
        end
    end
end