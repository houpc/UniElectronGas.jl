using ElectronLiquid, UniElectronGas
using JLD2, DelimitedFiles

include("../input.jl")

if isLayered2D
    const filename = "./data_freeE_layered2d_sg.jld2"
    const savefilename = "freeE_$(dim)d_layered2d.dat"
    const _parafilename = "para_wn_1minus0_layered2d.csv"
else
    const filename = "./data$(dim)d_freeE.jld2"
    const savefilename = "freeE_$(dim)d.dat"
    const _parafilename = "para_wn_1minus0.csv"
end
# const filename = "data3d/data$(dim)d_freeE.jld2"
const filename_E0 = "E0_$(dim)d.txt"

function free_energy_0(para)
    if para.dim == 2
        return para.basic.n * para.EF / 2
    elseif para.dim == 3
        return para.basic.n * para.EF * 3 / 5
    else
        error("unknown dimension")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")
    paras = [UEG.paraid(ParaMC(k)) for k in keys(f)]

    # E0_data = readdlm(filename_E0)
    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, isFock=isFock, dim=dim)
        kF = para.kF
        if UEG.paraid(para) in paras
            println("working on ", UEG.paraid(para))
            energy0 = free_energy_0(para)
            energy = UniElectronGas.getEnergy(para, filename, [energy0, 0.0]; parafile=_parafilename)
            push!(results, Any[_rs, _beta, _mass2, _order, spin, energy...])
        else
            println("missing ", UEG.paraid(para), " in jld2 MC data file.")
        end
    end

    if isSave
        open(savefilename, "a+") do io
            writedlm(io, results)
        end
    end
end