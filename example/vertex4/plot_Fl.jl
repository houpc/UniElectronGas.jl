using PyPlot
using DelimitedFiles

dim = 3
spin = 2
# rs = [0.5, 1.0, 4.0]
# rs = [1.0]
rs = [3.0,]
ells = [0,]
# symmetry = true
symmetry = false
# mass2 = [1.0, 1.5, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
# mass2 = [1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
# mass2 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
mass2 = [1.0, 1.5, 2.0, 2.5, 3.0]
Fs = [-0.0,]
beta = [100.0]
order = [5,]
# const fileName = symmetry ? "Fsl_$(dim)d.dat" : "Fal_$(dim)d.dat"
const fileName = "gayuk_3d.dat"

cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function plot_convergence(ver4, errors, _mass2=mass2, maxOrder=order[1]; rs=rs[1], beta=beta[1])
    style = PyPlot.matplotlib."style"
    # style.use(["science", "std-colors"])
    style.use(["science"])
    # color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["magenta"], cdict["grey"], "black"]
    color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["grey"], "black"]
    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(6, 4))

    x = collect(1:maxOrder)
    for (idx, yval) in enumerate(ver4)
        yerr = errors[idx]
        errorbar(x, yval, yerr=yerr, color=color[idx], capsize=4, fmt="o", markerfacecolor="none", label="$(_mass2[idx])")
    end
    # xlim(0.8, 4.5)
    xlabel("Order")
    ylabel("\$m^*/m\$")
    legend(title="mass2")
    title("rs=$rs, beta=$beta")
    # title(r"$r_s=$(r_, beta=$")
    savefig("ver4$(dim)d_rs$(rs)_beta$(beta)_spin1_conv_fitted.pdf")
end

function plot_convergence_v1(ver4, errors, _mass2=mass2, maxOrder=order[1]; rs=rs[1], beta=beta[1], ell=ells[1])
    style = PyPlot.matplotlib."style"
    # style.use(["science", "std-colors"])
    style.use(["science"])
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
        yval, yerr = ver4[o], errors[o]
        # errorbar(x[2:end], yval[2:end], yerr=yerr[2:end], color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="$o")
        errorbar(x, yval, yerr=yerr, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="$o")
    end
    # xlim(0.8, 4.5)
    xlabel("lambda")
    ylabel("\$F^s_$(ell)\$")
    legend(title="order")
    title("rs=$rs, beta=$beta")
    # title(r"$r_s=$(r_, beta=$")
    if symmetry
        savefig("Fs$(dim)d_rs$(rs)_beta$(beta)_l$(ell)_conv.pdf")
    else
        savefig("Fa$(dim)d_rs$(rs)_beta$(beta)_l$(ell)_conv.pdf")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__

    ver4_data = readdlm(fileName)
    num_data = size(ver4_data)[1]

    for _l in ells
        ver4_total, error_total, mass2_total = [], [], []
        for (_rs, _beta, _mass2) in Iterators.product(rs, beta, mass2)
            idx = 0
            # println(_rs, _beta, _mass2)
            for i in 1:num_data
                if ver4_data[i, 1:3] == [_rs, _beta, _mass2] && ver4_data[i, 5] == _l
                    idx = i
                    break
                end
            end
            idx == 0 && continue
            _order = ver4_data[idx, 4]
            ver4, error = [], [], []
            for o in 1:_order
                push!(ver4, ver4_data[idx, 3o+3])
                push!(error, ver4_data[idx, 3o+5])
            end
            push!(ver4_total, ver4)
            push!(error_total, error)
            push!(mass2_total, ver4_data[idx, 3])
        end
        println(ver4_total)
        println(error_total)
        println(mass2_total)
        # plot_convergence(ver4_total, error_total, mass2_total)

        ver4_total = hcat(ver4_total...)
        error_total = hcat(error_total...)
        ver4_order, error_order = [], []
        for o in 1:order[1]
            push!(ver4_order, ver4_total[o, :])
            push!(error_order, error_total[o, :])
        end

        plot_convergence_v1(ver4_order, error_order, mass2_total, ell=_l)
    end
end