using ElectronLiquid, UniElectronGas
using JLD2, DelimitedFiles

dim = 3
rs = [1.0]
# mass2 = [0.5, 1.0, 2.0, 2.5, 3.0, 4.0]
# mass2 = [1.6, 1.8, 2.0, 2.2]
mass2 = [2.3]
Fs = [-0.0]
beta = [80.0]
# beta = [25.0, 40.0, 80.0]
order = [4]
isDynamic = false
isFock = false
spinPolarPara = 0.0

const parafilename = "para_wn_1minus0.csv"
const filename = "data3d_freeE.jld2"
# const filename_E0 = "E0_3d.dat"
const filename_E0 = "E0_3d.txt"
const savefilename = "freeE_$(dim)d.dat"

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")
    E0_data = readdlm(filename_E0)
    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, isFock=isFock, dim=dim)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                println(UEG.paraid(para))
                idx = 0
                for i in 1:size(E0_data)[1]
                    if E0_data[i, 1:2] == [_rs, _beta]
                        idx = i
                        break
                    end
                end
                if idx == 0
                    energy = UniElectronGas.getEnergy(para, filename)
                else
                    energy = UniElectronGas.getEnergy(para, filename, E0_data[idx, 4:5])
                end
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