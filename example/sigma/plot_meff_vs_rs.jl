using CurveFit
using DelimitedFiles
using Measurements
using PyCall
using PyPlot

@pyimport scienceplots  # `import scienceplots` is required as of version 2.0.0
@pyimport scipy.interpolate as interp

function spline(x, y, e; xmin=x[1], xmax=x[end])
    # generate knots with spline without constraints
    w = 1.0 ./ e
    spl = interp.UnivariateSpline(x, y; w=w, k=3)
    __x = collect(LinRange(xmin, xmax, 100))
    yfit = spl(__x)
    return __x, yfit
end

function spline_with_bc(x, y, e; xmin=0.0, xmax=x[end])
    _x, _y = deepcopy(x), deepcopy(y)
    _w = 1.0 ./ e

    #enforce left boundary condition: m*/m → 1 as rs → 0
    rescale = 10000
    pushfirst!(_x, 0.0)
    pushfirst!(_y, 1.0)
    pushfirst!(_w, _w[1] * rescale)

    # generate knots with spline without constraints
    spl = interp.UnivariateSpline(_x, _y; w=_w, k=3)
    __x = collect(LinRange(0.0, x[end], 1000))
    yfit = spl(__x)
    return __x, yfit
end

dim = 3
spin = 2
Fs = -0.0
beta = 40.0
maxOrder = 5
rslist = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

# meff_data = Dict(
#     1.0 => [0.93574 ± 0.0001, 0.94731 ± 0.00011, 0.94992 ± 0.00013, 0.95129 ± 0.00025, 0.95109 ± 0.00092],
#     2.0 => [0.97009 ± 0.00027, 0.95709 ± 0.00035, 0.95099 ± 0.00041, 0.9504 ± 0.00066, 0.9516 ± 0.0019],
#     3.0 => [0.94561 ± 0.00012, 0.95116 ± 0.00021, 0.95442 ± 0.00023, 0.95754 ± 0.00039, 0.96014 ± 0.00054],
#     4.0 => [0.95522 ± 0.00012, 0.96049 ± 0.00023, 0.96212 ± 0.00029, 0.96639 ± 0.00054, 0.968 ± 0.0014],
#     5.0 => [0.95845 ± 0.00029, 0.96925 ± 0.00068, 0.9696 ± 0.001, 0.9752 ± 0.0015, 0.9766 ± 0.0052],
# )

meff_data = Dict(
    # Fixed lambda for rs < 3
    0.5 => [0.94957 ± 0.00012, 0.95751 ± 0.00013, 0.95808 ± 0.00013, 0.95839 ± 0.00016, 0.95861 ± 0.00041],
    1.0 => [0.93574 ± 0.0001, 0.94731 ± 0.00011, 0.94992 ± 0.00013, 0.95129 ± 0.00025, 0.95109 ± 0.00092],
    2.0 => [0.97009 ± 0.00027, 0.95709 ± 0.00035, 0.95099 ± 0.00041, 0.9504 ± 0.00066, 0.9516 ± 0.0019],
    # Minimal sensitivity for rs ≥ 3
    # NOTE: N = 1 has a plateau at λ = ∞ => (m*/m)₁ = 1
    3.0 => [1.0 ± 0.0, 0.95116 ± 0.00021, 0.95442 ± 0.00023, 0.95754 ± 0.00039, 0.96014 ± 0.00054],
    4.0 => [1.0 ± 0.0, 0.96049 ± 0.00023, 0.96212 ± 0.00029, 0.96639 ± 0.00054, 0.968 ± 0.0014],
    5.0 => [1.0 ± 0.0, 0.96925 ± 0.00068, 0.9696 ± 0.001, 0.9752 ± 0.0015, 0.9766 ± 0.0052],
)
lambdas = Dict(
    0.5 => [3.5, 3.5, 3.5, 3.5, 3.5],
    1.0 => [1.75, 1.75, 1.75, 1.75, 1.75],
    2.0 => [2.0, 2.0, 2.0, 2.0, 2.0],
    3.0 => [0.75, 0.75, 1.0, 1.25, 1.75],
    4.0 => [0.625, 0.625, 0.75, 1.0, 1.125],
    5.0 => [0.5, 0.5, 0.625, 0.875, 0.875],
)

style = PyPlot.matplotlib."style"
style.use(["science", "std-colors"])
color = [
    "black",
    cdict["orange"],
    cdict["blue"],
    cdict["cyan"],
    cdict["magenta"],
    cdict["red"],
    # cdict["teal"],
]
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 16
rcParams["font.family"] = "Times New Roman"
figure(figsize=(6, 4))

for o in 1:maxOrder
    # for o in maxOrder:maxOrder
    yval = [meff_data[rs][o].val for rs in rslist]
    yerr = [meff_data[rs][o].err for rs in rslist]
    println(rslist)
    println(yval)
    println(yerr)

    # Plot order-by-order lambda optima for rs > 2 at d = 3
    errorbar(
        rslist,
        yval,
        yerr=yerr,
        # color="black",
        color=color[o],
        capsize=4,
        fmt="o",
        markerfacecolor="none",
        label="\$N = $o\$",
        zorder=10 * o,
    )
    if o > 1
        xfit, yfit = spline_with_bc(rslist, yval, yerr)
        # xfit, yfit = spline(rslist, yval, yerr)
        plot(xfit, yfit; color=color[o], linestyle="--")
        # plot(xfit, yfit; color="black", linestyle="--")

        # Estimate the turning point (local minimum) of the mass ratio to max order
        if o == maxOrder
            println("\nTurning point: rs = $(xfit[argmin(yfit)])")
            println("Effective mass ratio at turning point: $(minimum(yfit))")
        end
    else
        axhline(1.0; color=color[o], linestyle="--")
    end
end
text(
    3.5,
    0.94,
    "\$\\beta \\hspace{0.1em} \\epsilon_F = $(beta)\$";
    fontsize=16
)
legend(; loc=(0.1, 0.575), ncol=2)
# xticks(range(0.0, 5.0, step=0.5))
xlim(-0.2, 5.2)
ylim(0.93, 1.005)
xlabel("\$r_s\$")
ylabel("\$m^\\star / m\$")
if spin == 2
    savefig("meff$(dim)d_beta$(beta)_vs_rs.pdf")
else
    savefig("meff$(dim)d_beta$(beta)_spin$(spin)_vs_rs.pdf")
end

# Plot meff convergence at each rs
for rs in rslist
    figure(figsize=(6, 4))
    orders = 1:maxOrder
    yval = Measurements.value.(meff_data[rs])
    yerr = Measurements.uncertainty.(meff_data[rs])
    println(yval)
    println(yerr)

    # Plot order-by-order lambda optima for rs > 2 at d = 3
    errorbar(
        orders,
        yval,
        yerr=yerr,
        color=cdict["blue"],
        capsize=4,
        fmt="o-",
        markerfacecolor="none",
    )
    if rs > 2
        for o in orders
            lambdastr = o == 1 ? "\\infty" : "$(lambdas[rs][o])"
            annotate(
                "\$$lambdastr\$",
                (o, yval[o]),
                xytext=(o + 0.1, yval[o] - 0.0025),
                fontsize=12,
            )
        end
    end

    # xfit, yfit = spline(orders, yval, yerr)
    # plot(xfit, yfit; color=cdict["blue"], linestyle="--")
    lambdastr = rs < 3 ? ",\\, \\lambda^\\star = $(lambdas[rs][1])" : ""
    yloc = 0.99
    if rs == 0.5
        yloc=0.9505
        xlim(0.75, 5.25)
        ylim(0.9485, 0.9605)
    elseif rs == 1.0
        yloc=0.9375
        xlim(0.75, 5.25)
        ylim(0.934, 0.956)
    elseif rs == 2.0
        yloc = 0.9675
        xlim(0.75, 5.25)
        ylim(0.949, 0.971)
    elseif rs == 3.0
        yloc = 0.99
        xlim(0.75, 5.65)
        ylim(0.945, nothing)
    elseif rs == 4.0
        yloc = 0.992
        xlim(0.75, 5.65)
        ylim(0.955, nothing)
    elseif rs == 5.0
        yloc = 0.994
        xlim(0.75, 5.65)
        ylim(0.965, nothing)
    end
    text(
        rs < 3 ? 2.0 : 3.0,
        yloc,
        "\$r_s = $(rs),\\, \\beta \\hspace{0.1em} \\epsilon_F = $(beta)$lambdastr\$";
        fontsize=16
    )
    xlabel("Perturbation order \$N\$")
    ylabel("\$m^\\star / m\$")
    if spin == 2
        savefig("meff$(dim)d_rs=$(rs)_beta$(beta)_vs_N.pdf")
    else
        savefig("meff$(dim)d_rs=$(rs)_beta$(beta)_spin$(spin)_vs_N.pdf")
    end
end
