using UniElectronGas, ElectronLiquid
using JLD2

dim = 3
rs = [4.0,]
mass2 = [0.5,]
Fs = [-0.0,]
beta = [25.0,]
order = [5,]
isDynamic = false
# dSigma/dk = (Sigma[kF_label+idx_dk] - Sigma[kF_label-idx_dk]) / (kgrid[kF_label+idx_dk] - kgrid[kF_label-idx_dk])
idx_dk = 1

const parafilename = "para_wn_1minus0.csv"
const filename = "data_K.jld2"

if abspath(PROGRAM_FILE) == @__FILE__
    f = jldopen(filename, "r")

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                meff = getMeff(para, filename, idx_dk)
                println(UEG.paraid(para))
                println("meff: ", meff)
            end
        end
    end
end