using PyPlot, PyCall
using DelimitedFiles
using LsqFit

# @pyimport scienceplots  # `import scienceplots` is required as of version 2.0.0
@pyimport scipy.interpolate as interp

include("../input.jl")
# dim = 2
# spin = 2
# spin = 1
# rs = [0.5, 1.0, 4.0]
# mass2 = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0]
# Fs = [-0.0,]
# beta = [25.0]
# order = [5,]
const fileName = spin == 2 ? "zfactor_$(dim)d.dat" : "zfactor_$(dim)d_spin$spin.dat"

cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function spline(x, y, e; xmin=x[1], xmax=x[end])
    # generate knots with spline without constraints
    w = 1.0 ./ e
    spl = interp.UnivariateSpline(x, y; w=w, k=3)
    __x = collect(LinRange(xmin, xmax, 100))
    yfit = spl(__x)
    return __x, yfit
end

function plot_conv(zfactor, errors, _mass2=mass2, maxOrder=order[1]; rs=rs[1], beta=beta[1])
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
    for o in 1:maxOrder
        yval, yerr = zfactor[o], errors[o]
        println("Order $o: ", x, yval, yerr)
        errorbar(x, yval, yerr=yerr, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="\$N=\$ $o")

        # yfit = curve_fit(Polynomial, x, yval, 4)
        if o < 5
            xfit, yfit = spline(x, yval, yerr)
            # plot(xgrid, yfit.(xgrid), color=color[o])
            plot(xfit, yfit, color=color[o])
        end
    end
    mres = zeros(3)

    legend_loc = 4
    if spin == 1
        xlim(0.3, 5.2)
        ylim(0.85, 1.01)
    end
    # axhspan(mres[1], mres[3], color=cdict["cyan"], alpha=0.5)
    # plot([0, 6], [mres[2], mres[2]], color="black", ls="--")

    xlabel("\$\\lambda\$ (Ry)")
    ylabel("\$z\$ factor")
    legend(fontsize=12, loc=legend_loc)
    # title("rs=$rs, beta=$beta")
    title("\$r_s=\$ $rs")
    if spin == 2
        savefig("zfactor$(dim)d_rs$(rs)_beta$(beta)_conv.pdf")
    else
        savefig("zfactor$(dim)d_rs$(rs)_beta$(beta)_spin$(spin)_conv.pdf")
    end
end

function plot_Tscaling(zfactor, errors, _beta=beta, maxOrder=order[1]; rs=rs[1], mass2=mass2[1])
    style = PyPlot.matplotlib."style"
    style.use(["science", "std-colors"])
    # color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["magenta"], cdict["grey"], "black"]
    color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], "black", cdict["cyan"], cdict["grey"]]
    #cmap = get_cmap("Paired")
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 16
    rcParams["font.family"] = "Times New Roman"
    figure(figsize=(6, 4))

    x = 1.0 ./ _beta
    perm = sortperm(x)
    for o in 1:maxOrder
        yval, yerr = zfactor[o], errors[o]
        println("Order $o: ", x, yval, yerr)
        errorbar(x, yval, yerr=yerr, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="\$N=\$ $o")

        # yfit = curve_fit(Polynomial, x, yval, 1)
        # xgrid = [0.0, 0.5]
        xgrid = LinRange(0, 0.1, 200)
        if o > 1
            # @. model(x, p) = p[1] + p[2] * sqrt(x)
            # @. model(x, p) = p[1] + p[2] * x + p[3] * x^2
            @. model(x, p) = p[1] + p[2] * x
            # @. model(x, p) = p[1] * log(x + p[2])

            # p0 = [0.95, 0.02, -0.1]
            p0 = [0.95, 0.02]
            fit_param, fit_param_err = similar(p0), similar(p0)
            for iter in 0:1
                fit = curve_fit(model, x[perm[1:end-iter]], yval[perm[1:end-iter]], 1 ./ yerr[perm[1:end-iter]] .^ 2, p0)
                fit_param = coef(fit)
                fit_param_err = stderror(fit)
                p0 = fit_param
                println("$iter-th fitting: ", fit_param)
                println(fit_param_err)
            end
            plot(xgrid, model(xgrid, fit_param), color=color[o])

            z_T0 = model(0.0, fit_param)
            errorbar([0.0], [z_T0], yerr=[fit_param_err[1] * 2], color=color[o], capsize=4, fmt="o")
        else
            plot(xgrid, ones(length(xgrid)), color=color[o])
            errorbar([0.0], [1.0], yerr=[0.0], color=color[o], capsize=4, fmt="o")
        end
        # xfit, yfit = spline(_beta, yval, yerr)
        # plot(1.0 ./ xfit, yfit, color=color[o])
    end
    mres = zeros(3)

    legend_loc = 4
    # if spin == 1
    #     xlim(0.3, 5.2)
    #     ylim(0.85, 1.01)
    # end
    # axhspan(mres[1], mres[3], color=cdict["cyan"], alpha=0.5)
    # plot([0, 6], [mres[2], mres[2]], color="black", ls="--")

    # xyscale("log")
    # yscale("log")
    xlim(-0.01, 0.11)
    # ylim(0.92, 1.01)
    ylim(0.91, 1.01)
    xlabel("\$T\$")
    ylabel("\$z\$ factor")
    legend(fontsize=12, loc=legend_loc)
    # title("rs=$rs, beta=$beta")
    title("\$r_s=\$ $rs, $(dim)D UEG")
    if spin == 2
        # savefig("zfactor$(dim)d_rs$(rs)_lam$(mass2)_tscaling_n1.pdf")
        savefig("zfactor$(dim)d_rs$(rs)_lam$(mass2)_tscaling.pdf")
    else
        # savefig("zfactor$(dim)d_rs$(rs)_lam$(mass2)_spin$(spin)_tscaling.pdf")
        savefig("zfactor$(dim)d_rs$(rs)_lam$(mass2)_spin$(spin)_tscaling.pdf")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__

    zfactor_data = readdlm(fileName)
    num_data = size(zfactor_data)[1]

    for _rs in rs
        zfactor_total, error_total, beta_total, mass2_total = [], [], [], []
        # for (_rs, _beta, _mass2) in Iterators.product(rs, beta, mass2)
        for (_beta, _mass2) in Iterators.product(beta, mass2)
            idx = 0
            # println(_rs, _beta, _mass2)
            for i in 1:num_data
                if zfactor_data[i, 1:3] == [_rs, _beta, _mass2]
                    # if zfactor_data[i, 1] == _rs && zfactor_data[i, 3] == _mass2
                    idx = i
                    break
                end
            end
            idx == 0 && continue
            _order = zfactor_data[idx, 4]
            zfactor, error = [], [], []
            for o in 1:_order
                push!(zfactor, zfactor_data[idx, 3o+2])
                push!(error, zfactor_data[idx, 3o+4])
            end
            if _order < order[1]
                push!(zfactor, 0.0)
                push!(error, 0.0)
            end
            push!(zfactor_total, zfactor)
            push!(error_total, error)
            push!(beta_total, zfactor_data[idx, 2])
            push!(mass2_total, zfactor_data[idx, 3])
        end
        println(zfactor_total)
        println(error_total)
        println(beta_total)
        println(mass2_total)

        zfactor_total = hcat(zfactor_total...)
        error_total = hcat(error_total...)
        zfactor_order, error_order = [], []
        for o in 1:order[1]
            push!(zfactor_order, zfactor_total[o, :])
            push!(error_order, error_total[o, :])
        end

        # plot_conv(zfactor_order, error_order, mass2_total, rs=_rs)
        plot_Tscaling(zfactor_order, error_order, beta_total, rs=_rs)
    end

    # plot_Tscaling([[1.0, 1, 1, 1, 1, 1], [1.000672, 0.998957, 0.99566, 0.992089, 0.98934, 0.987293]],
    #     [[3e-12, 5e-12, 1e-11, 2e-11, 4e-11, 8e-11], [1.2e-5, 2.2e-5, 3.5e-5, 4.9e-5, 7.2e-5, 9.6e-5]],
    #     [2.5, 5.0, 10.0, 20, 40, 80], 2)
end