using Test
using ElectronLiquid
using FeynmanDiagram
using FiniteDifferences
using Lehmann
using Measurements
using ElectronLiquid
using ElectronLiquid.CompositeGrids
using ElectronLiquid.UEG

function PP_interaction_dynamic(n, para::ParaMC, kamp=para.kF, kamp2=para.kF; kawargs...)
    kF = para.kF
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 16, 0.001, 16)
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * x * kamp * kamp2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = UEG.KO_W(q, n, para)
    end
    Wp *= para.NFstar / 4  # additional minus sign because the interaction is exchanged
    return Interp.integrate1D(Wp, xgrid)
end

function PPE_interaction_dynamic(n, para::ParaMC, kamp=para.kF, kamp2=para.kF; kawargs...)
    kF = para.kF
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 16, 0.001, 16)
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * x * kamp * kamp2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = UEG.KO_W(q, n, para)
    end
    Wp *= para.NFstar / 4  # additional minus sign because the interaction is exchanged
    return Interp.integrate1D(Wp, xgrid)
end

function rpa_from_rs(rs)
    para = UEG.ParaMC(rs=rs, beta=1000, Fs=0.0, mass2=1e-6, isDynamic=true)
    uu = -2 * PP_interaction_dynamic(0, para) + 2 * PP_interaction_dynamic(1, para)
    ud = -PP_interaction_dynamic(0, para) * 2
    # println("uu=$uu, ud=$ud")
    # println("rs=$rs, ", (1 * PP_interaction_dynamic(0, para) + 3 * PP_interaction_dynamic(1, para)) / 4)
    println("rs=$rs, ", (2uu - ud) / 2)
    return (2uu - ud) / 2, (-ud) / 2
end

function shift_u(u, wc1, wc2)
    # shift u from wc1 to wc2
    return u ./ (1 .+ u .* log(wc1 / wc2))#, uerr ./ (1 .+ u .* log(wc1 / wc2)) .^ 2
end


N = 32
# rslist = [0.0 + (4 / N) * n for n in 1:N]
# rslist = [1.0, 2.0, 3.0, 4.0, 5.7]
# rslist = [0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.7, 6.4, 5 * 1.6]
# rslist = [1, 2, 4, 5, 6, 10, 3.5]
rslist = [2.07, 2.3, 2.41, 2.66, 3.01, 5.7]
uclist = 0.0 .* rslist
mulist = 0.0 .* rslist
for i in 1:length(rslist)
    rs = rslist[i]
    uclist[i], mulist[i] = rpa_from_rs(rs)
    # para = UEG.ParaMC(rs=rs, beta=100, Fs=0.0, mass2=1e-6, isDynamic=true)
    # uu = -2 * PP_interaction_dynamic(0, para) + 2 * PP_interaction_dynamic(1, para)
    # ud = -PP_interaction_dynamic(0, para) * 2
    # # println("uu=$uu, ud=$ud")
    # # println("rs=$rs, ", (1 * PP_interaction_dynamic(0, para) + 3 * PP_interaction_dynamic(1, para)) / 4)
    # println("rs=$rs, ", (2uu - ud) / 2)
    # uclist[i] = (2uu - ud) / 2
    # mulist[i] = (-ud) / 2
end
println(uclist)
println(mulist)
println([shift_u(u, 1, 0.1) for u in uclist])
