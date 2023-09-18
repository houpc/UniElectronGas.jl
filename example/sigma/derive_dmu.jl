using UniElectronGas, ElectronLiquid
using JLD2, DelimitedFiles

dim = 3
# dim = 2
spin = 2
rs = [1.0]
mass2 = [2.0,]
# mass2 = [0.5, 1.0, 1.5, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
Fs = [-0.0,]
beta = [25.0,]
order = [4,]
isDynamic = false
spinPolarPara = 0.0

# const parafilename = "para_wn_1minus0.csv"
const parafilename = "para.csv"
# const filename = "./data2d_Z_v0.jld2"
const savefilename = "dmu_$(dim)d.dat"

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, spin=spin)
        dmu = UniElectronGas.getdmu(para)
        push!(results, Any[_rs, _beta, _mass2, _order, spinPolarPara, dmu...])
    end
    if isSave
        open(savefilename, "a+") do io
            writedlm(io, results)
        end
    end
end