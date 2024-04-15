using UniElectronGas, ElectronLiquid
using JLD2, DelimitedFiles

include("../input.jl")

### dSigma/dk = (Sigma[kF_label+idx_dk] - Sigma[kF_label-idx_dk]) / (kgrid[kF_label+idx_dk] - kgrid[kF_label-idx_dk])

const filename = "./data$(dim)d/data$(dim)d_dSig1.jld2"
const savefilename = spin == 2 ? "dSigdk_$(dim)d.dat" : "dSigdk_$(dim)d_spin$spin.dat"

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

                ngrid, kgrid, rSw_k, iSw_k = UniElectronGas.get_dSigmadk(para, filename; parafile=parafilename)
                meff = UniElectronGas.getMeff(para, [rSw_k[i][1] for i in eachindex(rSw_k)]; parafile=parafilename)

                # println(rSw_k .* para.me / kF)
                # meff, fit_p = UniElectronGas.getMeff(para, rSw_k, kgrid)
                # meff = UniElectronGas.getMeff(para, rSw_k, kgrid)
                println("m* / m = ", meff)
                # push!(results, Any[_rs, _beta, _mass2, _order, meff...])
                # writedlm("fit.dat", fit_p)
            end
        end
    end
    if isSave
        open(savefilename, "a+") do io
            writedlm(io, results)
        end
    end
end
