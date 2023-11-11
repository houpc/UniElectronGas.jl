using PyPlot
using DelimitedFiles

dim = 3
spin = 2
# rs = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
# rs = [1.0, 2.0, 3.0, 4.0, 5.0]
rs = [2.0]
Fs = -[0.223, 0.380, 0.516, 0.639, 0.752] * 0.0
# rs = [1.0]
ells = [0,]
# symmetry = true
symmetry = false
# mass2 = [1.0, 1.5, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0]
# mass2 = [1e-3,]
# mass2 = [2.0, 2.5, 3.0, 3.5]
mass2 = [2.5, 3.0]
# Fs = [-0.0,]
beta = [100.0]
order = [5,]
# const fileName = symmetry ? "gs_$(dim)d.dat" : "ga_$(dim)d.dat"
# const fileName = symmetry ? "gsko_$(dim)d.dat" : "gako_$(dim)d.dat"
# const fileName = symmetry ? "guu_$(dim)d.dat" : "gud_$(dim)d.dat"
# fileName = "gayuk_3d.dat"
fileName = "gsingyuk3_3d_01.dat"

cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function plot_convergence_v1(ver4, errors, _rs=rs, maxOrder=order[1]; mass2=mass2[1], beta=beta[1], ell=ells[1])
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
    x = _rs
    for o in 1:maxOrder
        yval, yerr = ver4[o], errors[o]
        # yval = sum.(ver4[i] for i in 1:o)
        # yerr = sqrt.(sum.(errors[i] .^ 2 for i in 1:o) ./ o)
        # errorbar(x[2:end], yval[2:end], yerr=yerr[2:end], color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="$o")
        errorbar(x, -yval, yerr=yerr, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="$o")
    end
    # ylim(0.0, 0.4)
    xlim(0.0, 5.1)
    xlabel("\$r_s\$")
    ylabel("\$-U_0\$")
    legend(title="order")
    title("beta=$beta")
    # title(r"$r_s=$(r_, beta=$")
    if symmetry
        savefig("gs$(dim)d_beta$(beta)_l$(ell)_conv.pdf")
    else
        savefig("ga$(dim)d_beta$(beta)_l$(ell)_conv.pdf")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__

    ver4_data = readdlm(fileName)
    num_data = size(ver4_data)[1]

    for _l in ells
        ver4_total, error_total, rs_total = [], [], []
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
            push!(rs_total, ver4_data[idx, 1])
        end
        println(ver4_total)
        println(error_total)
        println(rs_total)
        # plot_convergence(ver4_total, error_total, rs_total)

        ver4_total = hcat(ver4_total...)
        error_total = hcat(error_total...)
        ver4_order, error_order = [], []
        for o in 1:order[1]
            push!(ver4_order, ver4_total[o, :])
            push!(error_order, error_total[o, :])
        end

        plot_convergence_v1(ver4_order, error_order, rs_total, ell=_l)
    end
end