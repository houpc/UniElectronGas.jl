using Colors
using JLD2
using PyCall
using PyPlot
using ElectronGas, ElectronLiquid, Parameters
using Lehmann, GreenFunc, CompositeGrids, CurveFit

@pyimport numpy as np   # for saving/loading numpy data
@pyimport scienceplots  # for style "science"
@pyimport scipy.interpolate as interp

# Vibrant qualitative colour scheme from https://personal.sron.nl/~pault/
const cdict = Dict([
    "orange" => "#EE7733",
    "blue" => "#0077BB",
    "cyan" => "#33BBEE",
    "magenta" => "#EE3377",
    "red" => "#CC3311",
    "teal" => "#009988",
    "grey" => "#BBBBBB",
]);
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

function spline(x, y, e; xmin=0.0, xmax=x[end])
    # generate knots with spline without constraints
    w = 1.0 ./ e
    spl = interp.UnivariateSpline(x, y; w=w, k=3)
    __x = collect(LinRange(xmin, xmax, 1000))
    yfit = spl(__x)
    return __x, yfit
end

function plot_mvsrs(rs, meff_data, color, label, ls="-"; rs_LDL=nothing, meff_LDL=nothing, zorder=nothing, ax=plt.gca(), lw=1.5)
    # Add data in the low-density limit to the fit, if provided
    if !isnothing(rs_LDL) && !isnothing(meff_LDL)
        rslist = unique([rs; rs_LDL])
        mdatalist = unique([meff_data; meff_LDL])
        P = sortperm(rslist)
        rs = rslist[P]
        meff_data = mdatalist[P]
    end

    # mfitfunc = interp.PchipInterpolator(rs, meff_data)
    mfitfunc = interp.Akima1DInterpolator(rs, meff_data)
    xgrid = np.arange(0, 6.2, 0.01)
    # ax.scatter(rs, meff_data; color=color, marker="o")
    if isnothing(zorder)
        handle, = ax.plot(xgrid, mfitfunc(xgrid), ls=ls, color=color, label=label, lw=lw)
    else
        handle, = ax.plot(xgrid, mfitfunc(xgrid), ls=ls, color=color, label=label, zorder=zorder, lw=lw)
    end
    yfit = np.ma.masked_invalid(mfitfunc(xgrid))
    # print(yfit)
    # print("Turning point: rs = ", xgrid[np.argmin(yfit)])
    # print("Effective mass ratio at turning point: ", np.min(yfit))
    return handle
end

function main()
    # System parameters
    dim = 3
    rs = 2.0
    δK = 5e-6
    beta = 1000.0
    param = Parameter.rydbergUnit(1.0 / beta, rs, dim)
    
    max_steps = 100
    int_type = :rpa
    if int_type == :rpa
        wstring = "W_0"
    elseif int_type == :ko_const_p
        wstring = "W_+"
    elseif int_type == :ko_const_pm
        wstring = "W_{\\pm}"
    end

    # Load imaginary-time self-energy data
    s = []
    si = []
    local n_max
    jldopen("sigma_GW_$(dim)d/sigma_$(int_type).jl", "r") do f
        for i in 0:max_steps
            try
                push!(s, f["Σ_$i"])
                push!(si, f["Σ_ins_$i"])
            catch
                n_max = i - 1
                break
            end
        end
    end

    # Get static self-energies
    kkfs = []
    s_iw0s = []
    si_iw0s = []
    for j in eachindex(s)
        s_iw = dlr_to_imfreq(to_dlr(s[j]))
        si_iw = dlr_to_imfreq(to_dlr(si[j]))
        idx_iw0 = locate(s_iw.mesh[1], 0)
        push!(kkfs, s[j].mesh[2] / param.kF)
        push!(s_iw0s, s_iw[idx_iw0, :])
        push!(si_iw0s, si_iw[idx_iw0, :])
    end

    # DLR cutoffs
    kF = param.kF
    Euv, rtol = 1000 * param.EF, 1e-11
    maxK, minK = 20kF, 1e-8kF
    Nk, order = 12, 8

    # G0W0
    meff_g0w = SelfEnergy.massratio(param, s[1], si[1], δK)[1]
    zfactor_g0w = SelfEnergy.zfactor(param, s[1]; ngrid=[-1, 0])[1]
    # s0, si0 = SelfEnergy.G0W0(
    #     param;
    #     Euv=Euv,
    #     rtol=rtol,
    #     Nk=Nk,
    #     maxK=maxK,
    #     minK=minK,
    #     order=order,
    #     int_type=int_type,
    #     # Fs=int_type == :ko_const_p ? -Fs : -0.0,   # NOTE: NEFT uses opposite sign convention for F!
    #     # Fa=int_type == :ko_const_pm ? -Fa : -0.0,  # NOTE: NEFT uses opposite sign convention for F!
    # )
    # meff_g0w = SelfEnergy.massratio(param, s0, si0, δK)[1]
    # zfactor_g0w = SelfEnergy.zfactor(param, s0; ngrid=[-1, 0])[1]

    # GW0
    meff_gw = SelfEnergy.massratio(param, s[end], si[end], δK)[1]
    zfactor_gw = SelfEnergy.zfactor(param, s[end]; ngrid=[-1, 0])[1]

    # Plot the real and imaginary parts of the static self-energy
    fig, ax = plt.subplots(; nrows=2, sharex=true)
    fig.set_figheight(6)
    fig.set_figwidth(6)
    colors = reverse(hex.((sequential_palette(0, length(s); s=100))))
    for j in eachindex(s)
        if j == length(s) - 1
            plot_mvsrs(kkfs[j], real(s_iw0s[j]), "limegreen", "\$N_\\text{iter} = $n_max\$"; zorder=1000, lw=2, ax=ax[1])
            plot_mvsrs(kkfs[j], imag(s_iw0s[j]), "limegreen", ""; zorder=1000, lw=2, ax=ax[2])
            # ax[1].plot(kkfs[j], real(s_iw0s[j]); color="limegreen", lw=2, label="\$N_\\text{iter} = $n_max\$", zorder=1000)
            # ax[2].plot(kkfs[j], imag(s_iw0s[j]); color="limegreen", lw=2, zorder=1000)
        else
            plot_mvsrs(kkfs[j], real(s_iw0s[j]), "#$(colors[j])", ""; ax=ax[1])
            plot_mvsrs(kkfs[j], imag(s_iw0s[j]), "#$(colors[j])", ""; ax=ax[2])
            # ax[1].plot(kkfs[j], real(s_iw0s[j]); color="#$(colors[j])")
            # ax[2].plot(kkfs[j], imag(s_iw0s[j]); color="#$(colors[j])")
        end
    end

    ax[1].annotate(
        "\$r_s = $(Int(rs))\$",
        xy=(0.775, 0.225),
        xycoords="axes fraction",
        ha="center",
        va="center",
    )
    ax[2].annotate(
        "\$\\left(\\frac{m^*}{m}\\right)_{G_0 $wstring} = $(round(meff_g0w; sigdigits=4))\$",
        xy=(0.75, 0.75),
        xycoords="axes fraction",
        ha="center",
        va="center",
    )
    ax[2].annotate(
        "\$\\left(\\frac{m^*}{m}\\right)_{G $wstring} = $(round(meff_gw; sigdigits=4))\$",
        xy=(0.757, 0.6),
        xycoords="axes fraction",
        ha="center",
        va="center",
    )
    ax[2].annotate(
        "\$Z_{G_0 $wstring} = $(round(zfactor_g0w; sigdigits=4))\$",
        xy=(0.783, 0.45),
        xycoords="axes fraction",
        ha="center",
        va="center",
    )
    ax[2].annotate(
        "\$Z_{G $wstring} = $(round(zfactor_gw; sigdigits=4))\$",
        xy=(0.79, 0.3),
        xycoords="axes fraction",
        ha="center",
        va="center",
    )

    ax[1].set_ylim(-0.35, 0.25)
    ax[1].set_yticks(-0.3:0.1:0.2)
    ax[1].set_xlim(0, 5)
    ax[2].set_xlim(0, 5)
    ax[1].legend(loc="best")
    # ax[2].legend(loc="best")
    ax[1].set_ylabel("\$\\text{Re}\\Sigma_c(k, i\\omega_0)\$")
    ax[2].set_ylabel("\$\\text{Im}\\Sigma_c(k, i\\omega_0)\$")
    ax[2].set_xlabel("\$k / k_F\$")
    plt.tight_layout()
    fig.savefig("static_gw_self_energy_$(int_type).pdf")
end

main()
