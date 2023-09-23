using PyPlot
using DelimitedFiles
using CurveFit

dim = 2
spin = 2
# spin = 1
# rs = [0.5, 1.0, 4.0]
rs = [1.0]
# rs = [2.0]
# mass2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
mass2 = [1.0, 1.5, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
# mass2 = [1.5, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
Fs = [-0.0,]
beta = [25.0]
# order = [4,]
order = [5,]
const fileName = spin == 2 ? "meff_$(dim)d.dat" : "meff_$(dim)d_spin$spin.dat"

cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function plot_convergence(meff, errors, _mass2=mass2, maxOrder=order[1]; rs=rs[1], beta=beta[1])
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], "black", cdict["cyan"], cdict["grey"]]
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
    # xlim(0.8, 4.5)
    xlabel("Order")
    ylabel("\$m^*/m\$")
    legend(title="mass2")
    # title("rs=$rs, beta=$beta")
    title("rs=$rs")
    savefig("meff$(dim)d_rs$(rs)_beta$(beta)_spin1_conv_fitted.pdf")
end

function plot_convergence_v1(meff, errors, _mass2=mass2, maxOrder=order[1]; rs=rs[1], beta=beta[1])
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    # color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["magenta"], cdict["grey"], "black"]
    color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], "black", cdict["cyan"], cdict["grey"]]
    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(6, 4))

    x = _mass2
    xgrid = LinRange(0.0, 5.0, 100)
    for o in 1:maxOrder
        yval, yerr = meff[o], errors[o]
        errorbar(x, yval, yerr=yerr, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="$o")

        yfit = curve_fit(Polynomial, x, yval, 5)
        o < 5 && plot(xgrid, yfit.(xgrid), color=color[o])
    end
    xlim(0.8, 5.2)    #rs=1
    ylim(0.945, 0.992)
    # ylim(0.88, 1.0)
    # xlim(0.8, 4.7)    #rs=0.5
    # ylim(0.94, 0.98)
    xlabel("lambda")
    ylabel("\$m^*/m\$")
    legend(title="order", loc=2)
    # title("rs=$rs, beta=$beta")
    title("rs=$rs")
    if spin == 2
        savefig("meff$(dim)d_rs$(rs)_beta$(beta)_conv1.pdf")
    else
        savefig("meff$(dim)d_rs$(rs)_beta$(beta)_spin$(spin)_conv1.pdf")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__

    meff_data = readdlm(fileName)
    num_data = size(meff_data)[1]
    # num_data = 8

    # mass2_all = meff_data[:, 3]
    # meff_order, error_order = [], []
    # for o in 1:4
    #     push!(meff_order, meff_data[:, 3o+2])
    #     push!(error_order, meff_data[:, 3o+4])
    # end

    meff_total, error_total, mass2_total = [], [], []
    for (_rs, _beta, _mass2) in Iterators.product(rs, beta, mass2)
        idx = 0
        # println(_rs, _beta, _mass2)
        for i in 1:num_data
            # if meff_data[i, 1:3] == [_rs, _beta, _mass2]
            if meff_data[i, 1] == _rs && meff_data[i, 3] == _mass2
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
        if _order < order[1]
            push!(meff, 0.0)
            push!(error, 0.0)
        end
        push!(meff_total, meff)
        push!(error_total, error)
        push!(mass2_total, meff_data[idx, 3])
    end
    println(meff_total)
    println(error_total)
    println(mass2_total)
    # plot_convergence(meff_total, error_total, mass2_total)

    meff_total = hcat(meff_total...)
    error_total = hcat(error_total...)
    meff_order, error_order = [], []
    for o in 1:order[1]
        push!(meff_order, meff_total[o, :])
        push!(error_order, error_total[o, :])
    end

    plot_convergence_v1(meff_order, error_order, mass2_total)
end