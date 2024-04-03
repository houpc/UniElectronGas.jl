using UniElectronGas
using PyCall
using PyPlot
using DelimitedFiles
using CurveFit
using Measurements

@pyimport scienceplots  # `import scienceplots` is required as of version 2.0.0
@pyimport scipy.interpolate as interp

include("../input.jl")
res_directory = "./data"

const inverse_dispersion_ratio_filename = joinpath(res_directory, "inverse_dispersion_ratio_$(dim)d.dat")
const zinv_filename = joinpath(res_directory, "inverse_zfactor_$(dim)d.dat")
const meff_filename = joinpath(res_directory, "meff_$(dim)d.dat")

cdict = Dict([
    "blue" => "#0077BB",
    "cyan" => "#33BBEE",
    "teal" => "#009988",
    "orange" => "#EE7733",
    "red" => "#CC3311",
    "magenta" => "#EE3377",
    "grey" => "#BBBBBB",
])
style = PyPlot.matplotlib."style"
style.use(["science", "std-colors"])
const color = [
    "black",
    cdict["orange"],
    cdict["blue"],
    cdict["cyan"],
    cdict["magenta"],
    cdict["red"],
    cdict["teal"],
]
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 16
rcParams["mathtext.fontset"] = "cm"
rcParams["font.family"] = "Times New Roman"

# rₛ ↦ λ*(rₛ), d = 3
const fixed_lambda_optima_3d = Dict(
    0.5 => 3.5,
    1.0 => 1.75,
    2.0 => 2.0,
    3.0 => 1.5,
    4.0 => 1.125,
    5.0 => 1.125,
    6.0 => missing,
)

# Two lambda points to plot for convergence tests
# Where possible, we take: (1) the chosen optimum λ*(rₛ) and (2) the largest calculated λ
# NOTE: At rs = 5, λ* itself is the largest calculated λ, so we compare with smaller λ = 0.875
const lambdas_meff_convergence_plot_3d = Dict(
    0.5 => [3.5, 5.0],
    1.0 => [1.75, 2.0, 3.0, 3.5],
    2.0 => [2.0, 2.5],
    3.0 => [1.5, 2.0],
    4.0 => [1.125, 1.5],
    5.0 => [0.875, 1.125],
    6.0 => [missing, missing],
)

### rs = 0.5 ###
# rs = [0.5]
# mass2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]

### rs = 1 ###
# rs = [1.0]
# mass2 = [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0]

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

function spline(x, y, e; xmin=x[1], xmax=x[end])
    # generate knots with spline without constraints
    w = 1.0 ./ e
    spl = interp.UnivariateSpline(x, y; w=w, k=3)
    __x = collect(LinRange(xmin, xmax, 100))
    yfit = spl(__x)
    return __x, yfit
end

# polarstr = ispolarized ? "_GV_spin_polarized" : ""
polarstr = ""

function spline(x, y, e; xmin=0.0, xmax=x[end])
    # generate knots with spline without constraints
    w = 1.0 ./ e
    spl = interp.UnivariateSpline(x, y; w=w, k=3)
    __x = collect(LinRange(xmin, xmax, 1000))
    yfit = spl(__x)
    return __x, yfit
end

function spline_with_bc(x, y, e; xmin=0.0, xmax=x[end])
    _x, _y = deepcopy(x), deepcopy(y)
    _w = 1.0 ./ e

    #enforce left boundary condition: zero derivative at 1/N → 0
    rescale = 10000
    pushfirst!(_x, 0.0)
    pushfirst!(_y, y[1] / _w[1])
    pushfirst!(_w, _w[1] * rescale)

    # generate knots with spline without constraints
    spl = interp.UnivariateSpline(_x, _y; w=_w, k=3)
    __x = collect(LinRange(xmin, xmax, 1000))
    yfit = spl(__x)
    return __x, yfit
end

# NOTE: assumes the following row format: 
#       | rs | beta | mass2 | order | mean_1 ± error_1 | ... | mean_N ± error_N |
function load_from_dlm(filename, mass2; rs=rs[1], beta=beta[1])
    data = readdlm(filename)
    num_data = size(data)[1]
    # num_data = 8

    idx = 0
    for i in 1:num_data
        if data[i, 1] == rs && data[i, 3] == mass2
            idx = i
            break
        end
    end
    @assert idx != 0 "Data for rs = $(rs), mass2 = $(mass2) not found in file $(filename)"
    _order = data[idx, 4]
    mean, error = [], [], []
    for o in 1:_order
        push!(mean, data[idx, 3o+2])
        push!(error, data[idx, 3o+4])
    end
    # # Add padding values if this mass2 run has N < maxOrder
    npad = max(0, order[1] - _order)
    append!(mean, repeat([missing], npad))
    append!(error, repeat([missing], npad))

    mean_total = mean
    error_total = error
    mass2_total = data[idx, 3]
    println(mean_total)
    println(error_total)
    println(mass2_total)
    return mean_total, error_total, mass2_total
end
function load_from_dlm(filename; rs=rs[1], beta=beta[1], sortby="order")
    @assert sortby in ["order", "mass2"]

    data = readdlm(filename)
    num_data = size(data)[1]
    # num_data = 8

    mean_total, error_total, mass2_total = [], [], []
    for _mass2 in mass2
        idx = 0
        for i in 1:num_data
            if data[i, 1] == rs && data[i, 3] == _mass2
                idx = i
                break
            end
        end
        idx == 0 && continue
        _order = data[idx, 4]
        mean, error = [], [], []
        for o in 1:_order
            push!(mean, data[idx, 3o+2])
            push!(error, data[idx, 3o+4])
        end
        # # Add padding values if this mass2 run has N < maxOrder
        npad = max(0, order[1] - _order)
        append!(mean, repeat([missing], npad))
        append!(error, repeat([missing], npad))
        # Add results for this mass2 to lists
        push!(mean_total, mean)
        push!(error_total, error)
        push!(mass2_total, data[idx, 3])
    end

    println(mean_total)
    println(error_total)
    println(mass2_total)

    if sortby == "order"
        mt = hcat(mean_total...)
        et = hcat(error_total...)
        mean_order, error_order = [], []
        for o in 1:order[1]
            push!(mean_order, mt[o, :])
            push!(error_order, et[o, :])
        end
        return mean_order, error_order, mass2_total
    else # sortby == "mass2"
        return mean_total, error_total, mass2_total
    end
end

function plot_all_order_convergence(; plot_rs=[1.0, 5.0])
    @assert all(rs ∈ keys(fixed_lambda_optima_3d) for rs in plot_rs) "λ* not defined for all inputs rs = $(rs)"
    plot_lambda = [fixed_lambda_optima_3d[rs] for rs in plot_rs]

    figure(figsize=(6, 4))
    # for (j, (rs, lambda)) in enumerate(zip(plot_rs, plot_lambda))
    for (rs, lambda) in zip(plot_rs, plot_lambda)
        j = 2
        filenames = [
            inverse_dispersion_ratio_filename,
            zinv_filename,
            # meff_filename,
            # dispersion_ratio_filename,
        ]
        labels = [
            "\$(1 + \\delta m)^{-1}\$",
            "\$(1 - \\delta s)\$",
            # "\$m^\\star / m\$",
            # "\$(1 + \\delta m)\$",
        ]
        colors = [
            [cdict["orange"], cdict["red"]],
            [cdict["cyan"], cdict["blue"]],
            ["gray", "black"],
        ]
        axhline(1.0; linestyle="dotted", color=cdict["grey"], linewidth=1, zorder=-10)
        for (i, (filename, label)) in enumerate(zip(filenames, labels))
            means, errors, lambda = load_from_dlm(filename, lambda; rs=rs)

            valid_means = collect(skipmissing(means))
            valid_errors = collect(skipmissing(errors))
            x = reverse(1.0 ./ collect(eachindex(valid_means)))
            yval = reverse(valid_means)
            yerr = reverse(valid_errors)
            errorbar(
                x,
                yval,
                yerr=yerr,
                color=colors[i][j],
                capsize=4,
                fmt="o--",
                label=label,
                zorder=10 * i,
            )
        end

        label = "\$m^\\star / m\$"
        means, errors, lambda = load_from_dlm(meff_filename, lambda; rs=rs)
        i = 3

        valid_means = collect(skipmissing(means))
        valid_errors = collect(skipmissing(errors))
        x = reverse(1.0 ./ collect(eachindex(valid_means)))
        yval = reverse(valid_means)
        yerr = reverse(valid_errors)
        errorbar(
            x,
            yval,
            yerr=yerr,
            color=colors[i][j],
            capsize=4,
            fmt="o",
            label=label,
            zorder=10 * i,
        )

        # N = ∞ error estimate: difference between N = 4 and N = 5 means plus N = 5 error
        extrapolated_error = valid_errors[end] + abs(valid_means[end] - valid_means[end-1])

        # Plot extrapolated value at N = ∞
        _x, _y = spline(x, yval, yerr)
        plot(_x, _y, color=colors[i][j], linestyle="--")
        errorbar([0], [0.955], yerr=[0.008],
            # _x[1],
            # _y[1],
            # yerr=extrapolated_error,
            capsize=4,
            fmt="o",
            color=colors[i][j],
            zorder=10 * i,
        )
    end
    xlabel("\$1/N\$")
    # ywin = 0.11
    # ylim(nothing, 1 + ywin)
    ylim(0.85, 1.11)
    # text(
    #     0.3,
    #     1 + 0.825 * ywin,
    #     "\$r_s = $(plot_rs[1])\$";
    #     fontsize=16
    # )
    # text(
    #     0.75,
    #     1 + 0.825 * ywin,
    #     "\$r_s = $(plot_rs[2])\$";
    #     fontsize=16
    # )
    # legend(; loc=(0.2, 0.525), ncol=2, columnspacing=0.45)
    legend(; loc=1, fontsize=14)
    savefig("meff$(dim)d_rs$(plot_rs[1])_beta$(beta[1])_spin$(spin)_fitted$(polarstr)_with_cancellations_vs_Ninv.pdf")
end

function plot_meff_order_convergence(maxOrder=order[1]; rs=rs[1], beta=beta[1])
    meff_means, meff_errors, _mass2 = load_from_dlm(meff_filename; sortby="mass2")
    println("meff_means")
    println(meff_means)

    # Plot convergence for lambda points with max order data
    figure(figsize=(6, 4))
    x = reverse(1 ./ collect(1:maxOrder))
    ic = 1

    for (i, lambda) in enumerate(_mass2)
        any(ismissing.(meff_means[i])) && continue
        lambda ∉ lambdas_meff_convergence_plot_3d[rs] && continue

        yval = reverse(meff_means[i])
        yerr = reverse(meff_errors[i])
        errorbar(
            x,
            yval,
            yerr=yerr,
            capsize=4,
            fmt="o",
            color=color[ic],
            markerfacecolor="none",
            label="\$\\lambda = $(_mass2[i])\$",
            zorder=10 * i,
        )

        # N = ∞ error estimate: difference between N = 4 and N = 5 means plus N = 5 error
        extrapolated_error = meff_errors[i][end] + abs(meff_means[i][end] - meff_means[i][end-1])

        # Plot extrapolated value at N = ∞
        _x, _y = spline(x, yval, yerr)
        plot(_x, _y, color=color[ic], linestyle="--")
        errorbar(
            _x[1],
            _y[1],
            yerr=extrapolated_error,
            capsize=4,
            fmt="o",
            color=color[ic],
            markerfacecolor="none",
            zorder=10 * i,
        )
        println("\n(λ = $(_mass2[i]))\n(m*/m)_∞ ≈ $(measurement(_y[1], extrapolated_error))")
        ic += 1

    end
    xlabel("\$1/N\$")
    ylabel("\$m^\\star / m\$")
    text(
        0.05,
        0.9725,
        "\$r_s = $(rs)\$";
        fontsize=16
    )
    ylim(0.932, 0.98)
    legend(; loc="upper right", ncol=2, columnspacing=0.45)
    savefig("meff$(dim)d_rs$(rs)_beta$(beta)_spin$(spin)_fitted$(polarstr)_vs_Ninv.pdf")
end

function plot_meff_lambda_convergence(maxOrder=order[1]; rs=rs[1], beta=beta[1])
    meff_means, meff_errors, _mass2 = load_from_dlm(meff_filename; sortby="order")

    figure(figsize=(6, 4))
    xmin_plot = Inf
    xmax_plot = -Inf
    xgrid = LinRange(0.0, 5.0, 100)
    for o in 1:maxOrder
        valid_meff = skipmissing(meff_means[o])
        valid_errors = skipmissing(meff_errors[o])

        idx_valid_mass2 = collect(eachindex(valid_meff))
        x = mass2[idx_valid_mass2]

        yval = collect(valid_meff)
        yerr = collect(valid_errors)

        if o == maxOrder
            println("\n(N = $o)\nλ = $x\nm*/m = $(measurement.(yval, yerr))\n")
        end

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

        xmin_plot = min(x[1], xmin_plot)
        xmax_plot = max(x[end], xmax_plot)
        if o < 5
            xfit, yfit = spline(x, yval, yerr)
            plot(xfit, yfit; color=color[o], linestyle="--")
        end
    end
    if dim == 3
        if rs < 1
            xpad = 0.2
        elseif rs < 3
            xpad = 0.1
        elseif rs < 6
            xpad = 0.05
        else
            xpad = 0.2
        end
        xlim(xmin_plot - xpad, xmax_plot + xpad)
        ylim(0.855, 1.0)
        # Plot fixed lambda optima for rs = 1, 2 at d = 3
        if dim == 3 && rs in keys(fixed_lambda_optima_3d)
            if ispolarized
                lambda_optimum = fixed_lambda_optima_3d_GV_spin_polarized[rs]
            else
                lambda_optimum = fixed_lambda_optima_3d[rs]
            end
            axvline(lambda_optimum; linestyle="-", color="dimgray", zorder=-10)
        end
        if rs == 0.5
            xloc = 4.0
            yloc = 0.9825
            ylim(0.865, 1.005)
        elseif rs == 1.0
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
            xloc = 0.55
            yloc = 1.0025
            ylim(0.915, 1.02)
        elseif rs == 6.0
            xloc = 2.5
            yloc = 0.985
            ylim(0.965, 1.005)
        end
        xmin, xmax = xlim()
        ymin, ymax = ylim()
        if rs < 1
            xstep = 1.0
        elseif rs < 4
            xstep = 0.5
        else
            xstep = 0.25
        end
        big_xticks = collect(range(0.0, 7.0, step=xstep))
        big_yticks = collect(range(0.8, 1.2, step=0.025))
    end
    text(
        xloc,
        yloc,
        "\$r_s = $(rs),\\, \\beta \\hspace{0.1em} \\epsilon_F = $(beta)\$";
        fontsize=16
    )
    if rs == 0.5
        columnspacing = 0.45
    elseif rs == 5.0
        columnspacing = 1.8
    else
        columnspacing = 0.9
    end
    legend(; loc="lower right", ncol=2, columnspacing=columnspacing)
    xlabel("\$\\lambda\$ (Ry)")
    ylabel("\$m^\\star / m\$")
    if spin == 2
        savefig("meff$(dim)d_rs$(rs)_beta$(beta)_fitted$(polarstr)_vs_lambda.pdf")
    else
        savefig("meff$(dim)d_rs$(rs)_beta$(beta)_spin$(spin)_fitted$(polarstr)_vs_lambda.pdf")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_all_order_convergence(plot_rs=[1.0])
    plot_meff_order_convergence()
    # plot_meff_lambda_convergence()
end