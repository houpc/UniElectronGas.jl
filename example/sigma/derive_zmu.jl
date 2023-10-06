using UniElectronGas, ElectronLiquid
using JLD2, DelimitedFiles

include("../input.jl")

if isLayered2D
    const parafilename = "para_wn_1minus0_layered2d.csv"
    const filename = "./data$(dim)d/data_Z_layered2d.jld2"
else
    const parafilename = "para_wn_1minus0.csv"
    const filename = "./data$(dim)d/data_Z.jld2"
end
# const filename = "./data2d_Z_v0.jld2"
# const filename = "./data$(dim)d/data$(dim)d_Z_beta80_rs$(rs[1]).jld2"
const savefilename = spin == 2 ? "zfactor_$(dim)d.dat" : "zfactor_$(dim)d_spin$spin.dat"


function zfactor_renorm(dz, dzinv; isRenorm=true)
    if isRenorm
        sumzinv = accumulate(+, dzinv)
        return @. 1.0 / (1.0 + sumzinv)
    else
        sumz = accumulate(+, dz)
        return @. 1.0 + sumz
    end
end

function process(para, datatuple, isSave)
    dz, dzinv, dmu = UniElectronGas.get_dzmu(para, datatuple; parafile=parafilename, verbose=1, isSave)

    z = zfactor_renorm(dz, dzinv)
    println("Zfactor: ", z)
    return z
end

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")
    println(filename)
    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, spin=spin)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                println("loading ... ", UEG.paraid(para))
                zfactor = process(para, f[key], isSave)
                push!(results, Any[_rs, _beta, _mass2, _order, zfactor...])
            end
        end
    end
    if isSave
        open(savefilename, "a+") do io
            writedlm(io, results)
        end
    end
end