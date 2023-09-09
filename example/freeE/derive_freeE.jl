using ElectronLiquid, UniElectronGas
using JLD2, DelimitedFiles

dim = 3
rs = [1.0]
# mass2 = [1.0, 2.0, 3.0]
mass2 = [0.6, 0.7, 0.8, 0.9]
Fs = [-0.0]
beta = [25.0]
order = [4]
isDynamic = false
isFock = false
spinPolarPara = 0.0

const parafilename = "para_wn_1minus0.csv"
const filename = "data_freeE.jld2"
const savefilename = "freeE_$(dim)d.dat"

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")
    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, isFock=isFock, dim=dim)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                println(UEG.paraid(para))
                energy = UniElectronGas.getEnergy(para, filename)
                push!(results, Any[_rs, _beta, _mass2, _order, spinPolarPara, energy...])
            end
        end
    end

    if isSave
        open(savefilename, "a+") do io
            writedlm(io, results)
        end
    end
end