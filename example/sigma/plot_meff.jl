using PyPlot
using DelimitedFiles

dim = 2
# rs = [0.5, 1.0, 4.0]
# rs = [0.5]
rs = [1.0]
# mass2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mass2 = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
Fs = [-0.0,]
# beta = [20.0, 25.0, 40.0, 80.0]
beta = [25.0]
order = [4,]
fileName = "meff_$(dim)d.dat"

cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function plot_convergence(meff, errors, _mass2=mass2, maxOrder=order[1]; rs=rs[1], beta=beta[1])
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    # color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["magenta"], cdict["grey"], "black"]
    color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["grey"], "black"]
    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(6, 4))

    x = collect(1:maxOrder)
    for (idx, yval) in enumerate(meff)
        yerr = errors[idx]
        errorbar(x, yval, yerr=yerr, color=color[idx], capsize=4, fmt="o", markerfacecolor="none", label="$(_mass2[idx])")
    end
    xlim(0.8, 4.5)
    xlabel("Order")
    ylabel("\$m^*/m\$")
    legend(title="mass2")
    title("rs=$rs, beta=$beta")
    # title(r"$r_s=$(r_, beta=$")
    savefig("meff$(dim)d_rs$(rs)_beta$(beta)_conv_fitted.pdf")
end

if abspath(PROGRAM_FILE) == @__FILE__

    meff_data = readdlm(fileName)
    num_data = size(meff_data)[1]
    # num_data = 8

    meff_total, error_total, mass2_total = [], [], []
    for (_rs, _beta, _mass2) in Iterators.product(rs, beta, mass2)
        idx = 0
        # println(_rs, _beta, _mass2)
        for i in 1:num_data
            if meff_data[i, 1:3] == [_rs, _beta, _mass2]
                idx = i
                break
            end
        end
        idx == 0 && continue
        _order = meff_data[idx, 4]
        meff, error = [], [], []
        for o in 1:_order
            push!(meff, meff_data[idx, 3o+2])
            push!(error, meff_data[idx, 3o+4])
        end
        push!(meff_total, meff)
        push!(error_total, error)
        push!(mass2_total, meff_data[idx, 3])
    end
    println(meff_total)
    println(error_total)
    println(mass2_total)
    plot_convergence(meff_total, error_total, mass2_total)
end