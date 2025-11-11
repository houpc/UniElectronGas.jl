using UniElectronGas
using PyCall
using PyPlot
using DelimitedFiles
using CurveFit
using Measurements

@pyimport scienceplots  # `import scienceplots` is required as of version 2.0.0
@pyimport scipy.interpolate as interp
@pyimport matplotlib.gridspec as gridspec

include("../input.jl")

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
    cdict["teal"],
    cdict["red"],
]
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 16
rcParams["mathtext.fontset"] = "cm"
rcParams["font.family"] = "Times New Roman"

modestr = ""
if spin != 2
    modestr *= "_spin$(spin)"
end
if ispolarized
    modestr *= "_polarized"
end
if isLayered2D
    modestr *= "_layered"
end

# Fixed lambda optima, rₛ ↦ λ*(rₛ), d = 3
const fixed_lambda_optima_3d = Dict(
    # 0.5 => 3.5,
    1.0 => 5.0,
    2.0 => 3.0,
    3.0 => 1.875,
    4.0 => 1.625,
    5.0 => 1.25,
    # 6.0 => 1.0,
)

# Lambda points to include in the convergence plots vs N
const lambdas_meff_convergence_plot_3d = Dict(
    # 0.5 => [3.5, 5.0],
    # 1.0 => [1.75, 3.5, 5.0, 6.5, 8.0],
    1.0 => [3.5, 5.0, 6.5, 8.0],
    # 1.0 => [3.5, 5.0, 6.5],
    2.0 => [2.0, 2.5, 3.0, 3.5],
    # 3.0 => [1.5, 1.75, 1.875, 2.0, 2.25],
    3.0 => [1.5, 1.75, 1.875, 2.0],
    # 4.0 => [1.375, 1.5, 1.625, 1.75, 2.0],
    4.0 => [1.375, 1.5, 1.625, 1.75],
    # 4.0 => [1.25, 1.375, 1.5, 1.625, 1.75],
    # 5.0 => [1.125, 1.25, 1.375, 1.5, 1.75],
    5.0 => [1.125, 1.25, 1.375, 1.5],
    # 6.0 => [0.75, 1.0, 1.25],
)

# Lambda points to include in the error estimation for the final
# m*/m estimates and error bars in the convergence plots vs N
const lambdas_meff_convergence_plot_converged_3d = Dict(
    # 0.5 => [3.5, 5.0],
    1.0 => [3.5, 5.0, 6.5],
    2.0 => [2.0, 2.5, 3.0],
    # 1.0 => [3.5, 5.0, 6.5, 8.0],
    # 2.0 => [2.0, 2.5, 3.0, 3.5],
    3.0 => [1.5, 1.75, 1.875, 2.0],
    4.0 => [1.375, 1.5, 1.625, 1.75],
    5.0 => [1.125, 1.25, 1.375, 1.5],
)

# NOTE: assumes the following row format: 
#       | rs | beta | mass2 | order | mean_1 ± error_1 | ... | mean_N ± error_N |
function load_from_dlm(filename, mass2; rs=rs[1], beta=beta[1], verbose=false)
    data = readdlm(filename)
    num_data = size(data)[1]
    idx = 0
    currMaxOrder = 0
    for i in 1:num_data
        if data[i, 1] == rs && data[i, 3] == mass2
            idx = i
            if data[i, 4] > currMaxOrder > 0
                println("(lambda = $(mass2)) Promoting from order $(currMaxOrder) to $(data[i, 4])")
                currMaxOrder = data[i, 4]
            end
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
    if verbose
        println(mean_total)
        println(error_total)
        println(mass2_total)
    end
    return mean_total, error_total, mass2_total
end
function load_from_dlm(filename; mass2=mass2, rs=rs[1], beta=beta[1], sortby="order", verbose=false)
    @assert sortby in ["order", "mass2"]

    data = readdlm(filename)
    num_data = size(data)[1]
    mean_total, error_total, mass2_total = [], [], []
    for _mass2 in mass2
        idx = 0
        currMaxOrder = 0
        for i in 1:num_data
            if data[i, 1] == rs && data[i, 3] == _mass2
                idx = i
                if data[i, 4] > currMaxOrder > 0
                    println("(lambda = $(_mass2)) Promoting from order $(currMaxOrder) to $(data[i, 4])")
                    currMaxOrder = data[i, 4]
                end
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
    if verbose
        println(mean_total)
        println(error_total)
        println(mass2_total)
    end

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

"""
Plot convergence of VDMC estimates for m*/m wrt perturbation order N.
Returns a dictionary of the final m*/m estimates with conservative error bars
determined by comparing a few lambda values demonstrating rapid convergence.
"""
function plot_meff_convergence_with_estimates(;
    beta=beta[1],
    plot_rs=range(1.0, 6.0),
    plot_lambdas=[lambdas_meff_convergence_plot_3d[rs] for rs in plot_rs],
)
    meff_estimates = Any[]
    num_rs = length(plot_rs)
    figure(figsize=(4 * num_rs, 4))
    for (i, (rs, lambdas)) in enumerate(zip(plot_rs, plot_lambdas))
        subplot(1, num_rs, i)
        best_mean = missing
        best_color = missing
        for (j, lambda) in enumerate(lambdas)
            means, errors, lambda = load_from_dlm(meff_dk_filename, lambda; rs=rs)
            valid_means = collect(skipmissing(means))
            yval = valid_means
            if lambda == fixed_lambda_optima_3d[rs]
                best_mean = yval[end]
                best_color = color[j+1]
                break
            end
        end
        max_err_through_N4 = 0.0
        for (j, lambda) in enumerate(lambdas)
            means, errors, lambda = load_from_dlm(meff_dk_filename, lambda; rs=rs)
            valid_means = collect(skipmissing(means))
            valid_errors = collect(skipmissing(errors))
            x = collect(eachindex(valid_means))
            yval = valid_means
            yerr = valid_errors
            # starstr = j == 1 ? "^*" : ""
            starstr = lambda == fixed_lambda_optima_3d[rs] ? "^*" : ""
            errorbar(
                x,
                yval,
                yerr=yerr,
                color=color[j+1],
                capsize=4,
                fmt="o--",
                # label="\$\\lambda = $lambda\$",
                label="\$\\lambda$starstr = $lambda\$",
                zorder=10 * j,
            )
            if lambda in lambdas_meff_convergence_plot_converged_3d[rs]
                d1 = abs(best_mean - yval[end])
                d2 = abs(best_mean - yval[end-1])
                d3 = abs(best_mean - yval[end-2])
                error_estimate = yerr[end] + max(d1, d2, d3)
                max_err_through_N4 = max(max_err_through_N4, error_estimate)
            end
        end
        # Rough estimate of total error using the last 3 orders
        # of a few lambda values demonstrating rapid convergence
        meff_estimate = measurement(best_mean, 1.1 * max_err_through_N4)
        push!(meff_estimates, Any[rs, meff_estimate])
        axhspan(
            best_mean - 1.1 * max_err_through_N4,
            best_mean + 1.1 * max_err_through_N4;
            color=best_color,
            alpha=0.2,
        )
        axhline(best_mean; linestyle="-", color=best_color, alpha=0.6)
        xticks(collect(1:order[1]))
        xlabel("Perturbation order \$N\$")
        if i == 1
            ylabel("\$m^* / m\$")
        end
        legend(; title="\$r_s = $(rs)\$", loc="best", fontsize=14, title_fontsize=14)
    end
    tight_layout()
    savefig("meff$(dim)d_beta$(beta)$(modestr)_vs_N.pdf")
    return meff_estimates
end

if abspath(PROGRAM_FILE) == @__FILE__
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    meff_estimates = plot_meff_convergence_with_estimates(plot_rs=[1.0, 2.0, 3.0, 4.0, 5.0])
    println("\nVDMC estimates for m*/m:")
    for (rs, meff) in sort(collect(meff_estimates))
        println("rs = $rs:\tm/m* ≈ $meff")
    end
    if isSave
        open(meff_estimates_filename, "a+") do io
            writedlm(io, meff_estimates)
        end
    end
end
