using DftFunctionals
using ElectronGas
using Measurements

function vxc_rs(rs_list)
    paras = [Parameter.atomicUnit(0.0001, rs) for rs in rs_list]
    ρ = reshape([p.n for p in paras], 1, :)
    ex = DftFunctional(:lda_x) # exchange
    co = DftFunctional(:lda_c_pw) # correlation, Perdew-Wang
    vex = potential_terms(ex, ρ).Vρ[1, :]
    vco = potential_terms(co, ρ).Vρ[1, :]
    vxc = vex .+ vco
    return vxc
end

function vc_rs(rs_list)
    paras = [Parameter.atomicUnit(0.001, rs) for rs in rs_list]
    ρ = reshape([p.n for p in paras], 1, :)
    co = DftFunctional(:lda_c_pw) # correlation, Perdew-Wang
    vco = potential_terms(co, ρ).Vρ[1, :]
    return vco
end

function muc(rs_list)
    paras = [Parameter.atomicUnit(0.001, rs) for rs in rs_list]
    mu = 0.0 .* rs_list
    for i in 1:length(rs_list)
        rs = rs_list[i]
        ρ = reshape([paras[i].n, paras[i].n * 1.01], 1, :)
        co = DftFunctional(:lda_c_pw) # correlation, Perdew-Wang
        eco = potential_terms(co, ρ).e
        # mu[i] = (eco[2] * ρ[2] - eco[1] * ρ[1]) / (ρ[2] - ρ[1])
        mu[i] = (eco[2] - eco[1]) / (ρ[2] - ρ[1])
    end
    return mu
end

function muxc(rs_list)
    paras = [Parameter.atomicUnit(0.001, rs) for rs in rs_list]
    mu = 0.0 .* rs_list
    for i in 1:length(rs_list)
        rs = rs_list[i]
        ρ = reshape([paras[i].n, paras[i].n * 0.001], 1, :)
        co = DftFunctional(:lda_c_pw) # correlation, Perdew-Wang
        eco = potential_terms(co, ρ).e
        ex = DftFunctional(:lda_x) # correlation, Perdew-Wang
        eex = potential_terms(ex, ρ).e
        exc = eco .+ eex
        # mu[i] = (eco[2] * ρ[2] - eco[1] * ρ[1]) / (ρ[2] - ρ[1])
        mu[i] = (exc[2] - exc[1]) / (ρ[2] - ρ[1])
    end
    return mu
end

rs_list = [1, 2, 4, 5, 10]

paras = [Parameter.atomicUnit(0.001, rs) for rs in rs_list]
ρ = reshape([p.n for p in paras], 1, :)
efs = [p.EF for p in paras]
println(efs)

ex = DftFunctional(:lda_x) # exchange
co = DftFunctional(:lda_c_pw) # correlation, Perdew-Wang
vex = potential_terms(ex, ρ).Vρ[1, :]
vco = potential_terms(co, ρ).Vρ[1, :]
eco = potential_terms(co, ρ).e
println("eco", eco)
vxc = vex .+ vco
muco = muc(rs_list)
println("muco", muco)
muxco = muxc(rs_list)
println("muxc", muxco)
println(vxc)
mu = (vxc ./ efs) .+ 1 # with EF included, divided by EF
println(mu)

# Holtzmann 2023 Tabel.S.2 SJ-GC-TABC-66
muqmc = [0.63 ± 0.01, 0.23 ± 0.01, -0.63 ± 0.01, -1.08 ± 0.01, -3.5 ± 0.02]

vxcqmc = (muqmc .- 1) .* efs
vcqmc = (muqmc .- 1) .* efs .- vex
println(vco)
println(vcqmc)
println(eco)

using PyPlot, LaTeXStrings
cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"])
style = PyPlot.matplotlib."style"
style.use(["science"])
color = [cdict["blue"], cdict["red"], cdict["teal"], cdict["orange"], cdict["cyan"], cdict["grey"], "black"]
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 16
rcParams["font.family"] = "Times New Roman"
figure(figsize=(6, 4))
xlabel(L"r_s")
# ylabel(L"$V$ (a.u.)")
ylabel(L"$V_{c}$ (a.u.)")

# plot(rs_list, vxc, label=L"PW")
rs_dense = [0.85 + 0.1 * n for n in 1:91]
# vxc_dense = vxc_rs(rs_dense)
vc_dense = vc_rs(rs_dense)
muc_dense = muc(rs_dense)

plot(rs_dense, vc_dense, label=L"PW, $V_c$")
plot(rs_dense, muc_dense, label=L"PW, $\mu_c$")
errorbar(rs_list, Measurements.value.(vcqmc), yerr=Measurements.uncertainty.(vcqmc), capsize=4,
    fmt="o", markerfacecolor="none", label=L"QMC, $\mu_c$")
# plot(rs_list, vxc, label=L"$V_{xc}$ PW")
# errorbar(rs_list, Measurements.value.(vxcqmc), yerr=Measurements.uncertainty.(vxcqmc), capsize=4,
#     fmt="o", markerfacecolor="none", label=L"$V_{xc}$ QMC")
# plot(rs_list, vco, label=L"$V_{c}$ PW")
# errorbar(rs_list, Measurements.value.(vcqmc), yerr=Measurements.uncertainty.(vcqmc), capsize=4,
#     fmt="o-", markerfacecolor="none", label=L"$V_{c}$ QMC")
legend()
xlim(0.0, 10.1)
# ylim(-0.1, 0.0)
# savefig("vc_from_mu.pdf")
savefig("vxc_from_mu.pdf")
# savefig("v_from_mu.pdf")