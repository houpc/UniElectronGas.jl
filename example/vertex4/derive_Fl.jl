using ElectronLiquid, UniElectronGas
using JLD2, DelimitedFiles

dim = 2
spin = 2
rs = [1.0]
mass2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
Fs = [-0.0]
beta = [25.0]
order = [3]
Nl = 2
isDynamic = false
isFock = false

const parafilename = "para_wn_1minus0.csv"
const filename = "data_ver4PH.jld2"
const savefilename1 = "Fsl_$(dim)d.dat"
const savefilename2 = "Fal_$(dim)d.dat"

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f = jldopen(filename, "r")
    results_s, results_a = Any[], Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, isFock=isFock, dim=dim, spin=spin)
        kF = para.kF
        for key in keys(f)
            loadpara = ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                println(UEG.paraid(para))
                data_Fs, data_Fa = UniElectronGas.getVer4PHl(para, filename)
                res_s = Any[_rs, _beta, _mass2, _order]
                for il in 1:Nl
                    push!(results_s, append!(Any[_rs, _beta, _mass2, _order, il-1], [real(data_Fs[o][il, 1]) for o in 1:_order]))
                    push!(results_a, append!(Any[_rs, _beta, _mass2, _order, il-1], [real(data_Fa[o][il, 1]) for o in 1:_order]))
                end
            end
        end
    end

    if isSave
        open(savefilename1, "a+") do io
            writedlm(io, results_s)
        end
        open(savefilename2, "a+") do io
            writedlm(io, results_a)
        end
    end
end