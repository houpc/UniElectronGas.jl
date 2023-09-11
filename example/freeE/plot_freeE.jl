using PyPlot
using DelimitedFiles

dim = 3
spin = 2
# rs = [0.5, 1.0, 4.0]
rs = [1.0]
# symmetry = true
symmetry = false
mass2 = [0.5, 1.0, 1.5, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
Fs = [-0.0,]
beta = [25.0]
order = [3,]
const fileName = "freeE_$(dim)d.dat"
const val_benchmark =
    cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function plot_convergence(freeE, errors, _mass2=mass2, maxOrder=order[1]; rs=rs[1], beta=beta[1])
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    # color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["magenta"], cdict["grey"], "black"]
    color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["grey"], "black"]
    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(6, 4))

    # x = collect(1:maxOrder)
    x = _mass2
    for o in 1:maxOrder
        yval, yerr = freeE[o], errors[o]
        # errorbar(x[2:end], yval[2:end], yerr=yerr[2:end], color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="$o")
        errorbar(x, yval, yerr=yerr, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="$o")
    end
    plot([0.2, 4.5], [1.174, 1.174])
    # xlim(0.8, 4.5)
    xlabel("lambda")
    ylabel("free energy")
    legend(title="order")
    title("rs=$rs, beta=$beta")
    # title(r"$r_s=$(r_, beta=$")
    savefig("freeE$(dim)d_rs$(rs)_beta$(beta)_conv.pdf")
end

if abspath(PROGRAM_FILE) == @__FILE__

    freeE_data = readdlm(fileName)
    num_data = size(freeE_data)[1]

    freeE_total, error_total, mass2_total = [], [], []
    for (_rs, _beta, _mass2) in Iterators.product(rs, beta, mass2)
        idx = 0
        # println(_rs, _beta, _mass2)
        for i in 1:num_data
            if freeE_data[i, 1:3] == [_rs, _beta, _mass2]
                idx = i
                break
            end
        end
        idx == 0 && continue
        _order = freeE_data[idx, 4]
        freeE, error = [], [], []
        for o in 1:_order
            push!(freeE, freeE_data[idx, 3o+6])
            push!(error, freeE_data[idx, 3o+8])
        end
        push!(freeE_total, freeE)
        push!(error_total, error)
        push!(mass2_total, freeE_data[idx, 3])
    end
    println(freeE_total)
    println(error_total)
    println(mass2_total)
    # plot_convergence(freeE_total, error_total, mass2_total)

    freeE_total = hcat(freeE_total...)
    error_total = hcat(error_total...)
    freeE_order, error_order = [], []
    for o in 1:order[1]
        push!(freeE_order, freeE_total[o, :])
        push!(error_order, error_total[o, :])
    end

    plot_convergence(freeE_order, error_order, mass2_total)
end