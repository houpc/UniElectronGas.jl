using PyCall
using PyPlot
using DelimitedFiles
using CurveFit
using Measurements

@pyimport scienceplots  # `import scienceplots` is required as of version 2.0.0
@pyimport scipy.interpolate as interp

# rₛ ↦ λ*(rₛ), d = 3
const fixed_lambda_optima_3d = Dict(
    1.0 => 1.75,
    2.0 => 2.0,
)
# rₛ ↦ λ*(rₛ, N), d = 3
const lambda_optima_3d = Dict(
    3.0 => [0.75, 0.75, 1.0, 1.25, 1.75],
    4.0 => [0.625, 0.625, 0.75, 1.0, 1.125],
)

lambdas_3d = Dict(
    1.0 => [1.75, 1.75, 1.75, 1.75, 1.75],
    2.0 => [2.0, 2.0, 2.0, 2.0, 2.0],
    3.0 => [0.75, 0.75, 1.0, 1.25, 1.75],
    4.0 => [0.625, 0.625, 0.75, 1.0, 1.125],
    5.0 => [0.5, 0.5, 0.625, 0.875, 0.875],
)

# rₛ ↦ λ*(rₛ), d = 3
const fixed_lambda_optima_3d_GV_spin_polarized = Dict(
    1.0 => 1.75,
)
lambdas_3d_GV_spin_polarized = Dict(
    1.0 => [1.75, 1.75, 1.75, 1.75, 1.75],
)

function spline(x, y, e; xmin=x[1], xmax=x[end])
    # generate knots with spline without constraints
    w = 1.0 ./ e
    spl = interp.UnivariateSpline(x, y; w=w, k=3)
    __x = collect(LinRange(xmin, xmax, 100))
    yfit = spl(__x)
    return __x, yfit
end

### rs = 1 ###
rs = [1.0]
mass2 = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0]

### rs = 2 ###
# rs = [2.0]
# mass2 = [0.5, 0.75, 1.0, 1.25, 1.5, 1.625, 1.75, 1.875, 2.0, 2.5, 3.0]

### rs = 3 ###
# rs = [3.0]
# mass2 = [0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.5, 1.75, 2.0]

### rs = 4 ###
# rs = [4.0]
# mass2 = [0.375, 0.5, 0.625, 0.75, 0.875, 1.0, 1.125, 1.25, 1.5, 2.0]

### rs = 5 ###
# rs = [5.0]
# mass2 = [0.375, 0.5, 0.625, 0.75, 0.8125, 0.875, 0.9375, 1.0, 1.125, 1.25, 1.5]

dim = 3
spin = 2
Fs = [-0.0,]
beta = [40.0]
order = [5,]

polarstr = ispolarized ? "_GV_spin_polarized" : ""

cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function spline(x, y, e; xmin=x[1], xmax=x[end])
    # generate knots with spline without constraints
    w = 1.0 ./ e
    spl = interp.UnivariateSpline(x, y; w=w, k=3)
    __x = collect(LinRange(xmin, xmax, 100))
    yfit = spl(__x)
    return __x, yfit
end

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
        errorbar(x, yval, yerr=yerr, color=color[idx], capsize=4, fmt="o", markerfacecolor="none", label="\$N=\$ $(_mass2[idx])")
    end
    xlabel("Order")
    ylabel("\$m^*/m\$")
    legend(title="mass2")
    # title("rs=$rs, beta=$beta")
    title("rs=$rs")
    savefig("meff$(dim)d_rs$(rs)_beta$(beta)_spin1_conv_fitted$(polarstr).pdf")
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

    xmin_plot = Inf
    xmax_plot = -Inf
    xgrid = LinRange(0.0, 5.0, 100)
    for o in 1:maxOrder
        valid_meff = skipmissing(meff[o])
        valid_errors = skipmissing(errors[o])

        idx_valid_mass2 = collect(eachindex(valid_meff))
        x = mass2[idx_valid_mass2]

        yval = collect(valid_meff)
        yerr = collect(valid_errors)

        # println(x)
        # println(yval)

        # yopt = yval[x .≈ lambdas_3d[rs][o]][1]
        # yerropt = yerr[x .≈ lambdas_3d[rs][o]][1]
        # println("N = $o, lambda = $(lambdas_3d[rs][o]):\tm*/m = $(yopt) ± $(yerropt)")

        # capsize = o == 5 ? 4 : nothing
        # errorbar(x, yval, yerr=yerr, color=color[o], markerfacecolor="none", capsize=capsize, fmt="o", label="\$N=$o\$", zorder=10 * o)

        # Plot order-by-order lambda optima for rs > 2 at d = 3
        errorbar(
            x,
            yval,
            yerr=yerr,
            color=color[o],
            capsize=4,
            fmt="o",
            markerfacecolor="none",
            label="\$N = $o\$",
            zorder=10 * o,
        )
        # if dim == 3 && rs > 2.0 && o > 1
        #     lambda_star_o = lambda_optima_3d[rs][o]
        #     errorbar(
        #         [lambda_star_o],
        #         yval[x.==lambda_star_o],
        #         yerr=yerr[x.==lambda_star_o],
        #         color=color[o],
        #         capsize=4,
        #         markersize=9,
        #         marker="*",
        #         # markerfacecolor="none",
        #         zorder=100 * o,
        #     )
        #     errorbar(
        #         x[x.!=lambda_star_o],
        #         yval[x.!=lambda_star_o],
        #         yerr=yerr[x.!=lambda_star_o],
        #         color=color[o],
        #         capsize=4,
        #         markersize=5,
        #         fmt="o",
        #         markerfacecolor="none",
        #         label="\$N = $o\$",
        #         zorder=10 * o,
        #     )
        # else
        #     errorbar(
        #         x,
        #         yval,
        #         yerr=yerr,
        #         color=color[o],
        #         capsize=4,
        #         markersize=5,
        #         fmt="o",
        #         markerfacecolor="none",
        #         label="\$N = $o\$",
        #         zorder=10 * o,
        #     )
        # end

        xmin_plot = min(x[1], xmin_plot)
        xmax_plot = max(x[end], xmax_plot)
        if o < 5
            xfit, yfit = spline(x, yval, yerr)
            plot(xfit, yfit; color=color[o], linestyle="--")
        end
    end
    if dim == 3
        xpad = rs < 3 ? 0.1 : 0.05
        xlim(xmin_plot - xpad, xmax_plot + xpad)
        ylim(0.855, 1.0)
        # Plot fixed lambda optima for rs = 1, 2 at d = 3
        if rs in [1.0, 2.0]
            if ispolarized
                lambda_optimum = fixed_lambda_optima_3d_GV_spin_polarized[rs]
            else
                lambda_optimum = fixed_lambda_optima_3d[rs]
            end
            axvline(lambda_optimum; linestyle="-", color="dimgray", zorder=-10)
        end
        if rs == 1.0
            if ispolarized
                xloc = 2.125
                yloc = 0.98
                ylim(0.83, 1.005)
                if spinPolarPara == 1.0
                    text(
                        0.8,
                        0.855,
                        "\$n_\\downarrow = 0\$";
                        fontsize=16
                    )
                else
                    text(
                        0.68,
                        0.98,
                        "\$\\frac{n_\\uparrow - n_\\downarrow}{n_\\uparrow + n_\\downarrow} = $spinPolarPara\$";
                        fontsize=16
                    )
                end
            else
                xloc = 2.125
                yloc = 0.98
                ylim(0.84, 1.005)
                text(
                    1.0,
                    0.855,
                    "\$n_\\uparrow = n_\\downarrow\$";
                    fontsize=16
                )
            end
        elseif rs == 2.0
            xloc = 0.75
            yloc = 0.9825
            ylim(0.865, 1.005)
        elseif rs == 3.0
            xloc = 0.6
            yloc = 0.8825
            ylim(0.865, 1.005)
        elseif rs == 4.0
            xloc = 0.625
            yloc = 0.905
            ylim(0.89, 1.005)
        elseif rs == 5.0
            xloc = 0.525
            yloc = 0.93
            ylim(0.915, 1.02)
        end
        xmin, xmax = xlim()
        ymin, ymax = ylim()
        xstep = rs < 4 ? 0.5 : 0.25
        big_xticks = collect(range(0.0, 5.0, step=xstep))
        big_yticks = collect(range(0.8, 1.2, step=0.025))
        xticks([t for t in big_xticks if xmin <= t <= xmax])
        yticks([t for t in big_yticks if ymin <= t <= ymax])
    end
    # xlim(0.8 * minimum(mass2), 1.2 * maximum(mass2))
    text(
        xloc,
        yloc,
        "\$r_s = $(rs),\\, \\beta \\hspace{0.1em} \\epsilon_F = $(beta)\$";
        fontsize=16
    )
    legend(; loc="lower right", ncol=2, columnspacing=0.9)
    xlabel("\$\\lambda\$ (Ry)")
    ylabel("\$m^\\star / m\$")
    # plt.tight_layout()
    # legend(title="order")
    # title("rs=$rs, beta=$beta")
    # title(r"$r_s=$(r_, beta=$")
    if spin == 2
        savefig("meff$(dim)d_rs$(rs)_beta$(beta)_fitted$(polarstr).pdf")
    else
        savefig("meff$(dim)d_rs$(rs)_beta$(beta)_spin$(spin)_fitted$(polarstr).pdf")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__

    meff_data = readdlm(meff_filename)
    num_data = size(meff_data)[1]
    # num_data = 8

    # mass2_all = meff_data[:, 3]
    # meff_order, error_order = [], []
    # for o in 1:4
    #     push!(meff_order, meff_data[:, 3o+2])
    #     push!(error_order, meff_data[:, 3o+4])
    # end

    for _rs in rs
        meff_total, error_total, mass2_total = [], [], []
        # for (_rs, _beta, _mass2) in Iterators.product(rs, beta, mass2)
        for (_beta, _mass2) in Iterators.product(beta, mass2)
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

        plot_convergence_v1(meff_order, error_order, mass2_total, rs=_rs)
    end
end