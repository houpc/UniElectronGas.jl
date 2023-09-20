using PyCall
using PyPlot
using DelimitedFiles

@pyimport scienceplots  # `import scienceplots` is required as of version 2.0.0

dim = 3
spin = 2
rs = [1.0]
mass2 = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0]
# mass2 = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0]
Fs = [-0.0,]
beta = [40.0]
order = [5,]
const fileName = spin == 2 ? "meff_$(dim)d.dat" : "meff_$(dim)d_spin$spin.dat"

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
    xlabel("Order")
    ylabel("\$m^*/m\$")
    legend(title="mass2")
    title("rs=$rs, beta=$beta")
    # title(r"$r_s=$(r_, beta=$")
    savefig("meff$(dim)d_rs$(rs)_beta$(beta)_spin1_conv_fitted.pdf")
end

function plot_convergence_v1(meff, errors, _mass2=mass2, maxOrder=order[1]; rs=rs[1], beta=beta[1])
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    # color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["magenta"], cdict["grey"], "black"]
    # color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["grey"], "black"]

    color = [
        "black",
        cdict["orange"],
        cdict["blue"],
        cdict["cyan"],
        cdict["magenta"],
        cdict["red"],
        # cdict["teal"],
    ]

    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["mathtext.fontset"] = "cm"
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(6, 4))

    for o in 1:maxOrder
        valid_meff = skipmissing(meff[o])
        valid_errors = skipmissing(errors[o])

        idx_valid_mass2 = collect(eachindex(valid_meff))
        x = mass2[idx_valid_mass2]

        yval = collect(valid_meff)
        yerr = collect(valid_errors)
        
        println(x)
        println(yval)

        # yval, yerr = meff[o], errors[o]

        # errorbar(x, yval, yerr=yerr, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="$o")
        capsize = o == 5 ? 4 : nothing
        errorbar(x, yval, yerr=yerr, color=color[o], capsize=capsize, markersize=3, fmt="o-", label="\$N=$o\$", zorder=10*o)
    end
    if dim == 3
        if rs == 1.0
            xloc = 2.0
            yloc = 0.971
            xlim(0.75, 4.0)
            ylim(0.865, 0.985)
            lambda_opt = 1.75
            axvline(lambda_opt; linestyle="--", color="dimgray", label="\$\\lambda^\\star = $(lambda_opt)\$")
        end
    end
    # xlim(0.8 * minimum(mass2), 1.2 * maximum(mass2))
    text(
        xloc,
        yloc,
        "\$r_s = $(rs),\\, \\beta \\hspace{0.1em} \\epsilon_F = $(beta)\$";
        fontsize=16
    )
    legend(; loc="lower right")
    xlabel("\$\\lambda\$ (Ry)")
    ylabel("\$m^\\star / m\$")
    # plt.tight_layout()
    # legend(title="order")
    # title("rs=$rs, beta=$beta")
    # title(r"$r_s=$(r_, beta=$")
    if spin == 2
        savefig("meff$(dim)d_rs$(rs)_beta$(beta)_fitted.pdf")
    else
        savefig("meff$(dim)d_rs$(rs)_beta$(beta)_spin$(spin)_fitted.pdf")
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
        # # Add padding values if this mass2 run has N < maxOrder
        npad = max(0, order[1] - _order)
        append!(meff, repeat([missing], npad))
        append!(error, repeat([missing], npad))
        # Add results for this mass2 to lists
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

    println(meff_total)
    println(error_total)

    meff_order, error_order = [], []
    for o in 1:order[1]
        push!(meff_order, meff_total[o, :])
        push!(error_order, error_total[o, :])
    end

    plot_convergence_v1(meff_order, error_order, mass2_total)
end