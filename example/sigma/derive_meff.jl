using UniElectronGas, ElectronLiquid
using JLD2, DelimitedFiles

# dim = 3
dim = 2
# rs = [0.5, 1.0, 4.0]
rs = [1.0]
mass2 = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
# mass2 = [3.0]
Fs = [-0.0,]
# beta = [20.0, 25.0, 40.0, 80.0]
beta = [25.0]
order = [4,]
isDynamic = false
### dSigma/dk = (Sigma[kF_label+idx_dk] - Sigma[kF_label-idx_dk]) / (kgrid[kF_label+idx_dk] - kgrid[kF_label-idx_dk])
# inds_dk = [1, 2, 3]

const parafilename = "para_wn_1minus0.csv"
# const filename = "./data2d/data$(dim)d_K.jld2"
# const filename = "./data$(dim)d_K.jld2"
const filename = "./data2d/data$(dim)d_K.jld2"

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")
    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                println(UEG.paraid(para))

                ngrid, kgrid, rSw_k, iSw_k = UniElectronGas.getSigma(para, filename)
                meff = UniElectronGas.getMeff(para, rSw_k, kgrid)
                println("m* / m = ", meff)
                # for idx_dk in inds_dk
                #     meff = UniElectronGas.getMeff(para, filename, idx_dk)
                #     println("dk index number: $idx_dk")
                #     println("m* / m = ", meff)
                #     println()
                #     push!(results, Any[_rs, _beta, _mass2, _order, idx_dk, meff...])
                # end
                push!(results, Any[_rs, _beta, _mass2, _order, meff...])
            end
        end
    end
    if isSave
        open("meff_$(dim)d.dat", "a+") do io
            writedlm(io, results)
        end
    end
end