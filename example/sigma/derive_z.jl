using UniElectronGas, ElectronLiquid
using JLD2, DelimitedFiles

include("../input.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(sigma_k_filename, "r")
    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, spin=spin)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                println(UEG.paraid(para))
                zinv = UniElectronGas.getZfactor(para; parafile=parafilename, isRenorm=false)
                println("1/z = ", zinv)
                push!(results, Any[_rs, _beta, _mass2, _order, zinv...])
            end
        end
    end
    if isSave
        open(zfactor_filename, "a+") do io
            writedlm(io, results)
        end
    end
end
