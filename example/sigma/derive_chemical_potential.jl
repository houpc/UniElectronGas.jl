using UniElectronGas, ElectronLiquid
using JLD2, DelimitedFiles

include("../input.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    results = Any[]
    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, spin=spin)
        try
            dmu = UniElectronGas.getdmu(para; parafile=parafilename)
            push!(results, Any[_rs, _beta, _mass2, _order, spinPolarPara, dmu...])
        catch
            println("No dmu data found for key $(UEG.short(para)) in CSV file $(parafilename)!")
        end
    end
    if isSave
        println(chemical_potential_filename)
        open(chemical_potential_filename, "a+") do io
            writedlm(io, results)
        end
    end
end
