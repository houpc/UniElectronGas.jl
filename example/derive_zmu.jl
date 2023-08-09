using UniElectronGas, ElectronLiquid
using JLD2

dim = 3
rs = [4.0,]
mass2 = [0.5,]
Fs = [-0.0,]
beta = [25.0,]
order = [5,]
isDynamic = false

const parafilename = "para_wn_1minus0.csv"
const filename = "data_Z.jld2"

function process(para, datatuple, isSave)
    if isSave
        UniElectronGas.save_zmu(para, datatuple; parafile=parafilename)
    end
    z = UniElectronGas.getZfactor(para, parafile=parafilename)
    println("Zfactor: ", z)
end

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                process(para, f[key], isSave)
            end
        end
    end
end