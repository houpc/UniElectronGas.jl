using ElectronLiquid, MCIntegration
using DelimitedFiles

spin = 2
dim = 2
rs = [1.0]
beta = [25.0]
neval = 1e6
isDynamic = false
spinPolarPara = spin / 2 - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1] 
filename = "E0_$(dim)d.txt"

function n(ω::T, β::T) where {T}
    return ω > T(0.0) ?
           exp(-ω * β) / (1 + exp(-ω * β)) :
           1.0 / (1 + exp(ω * β))
end

function integrand(x, config)
    para = config.userdata
    dim, β, me, μ = para.dim, para.β, para.me, para.μ
    k = x[1]

    ϵ = k^2 / (2me) - μ

    factor = dim == 2 ? 2π * k : 4π * k^2
    factor *= 1.0 / (2π)^dim

    return n(ϵ, β) * (ϵ + μ) * factor
end

for (_rs, _beta) in Iterators.product(rs, beta)
    para = UEG.ParaMC(rs=_rs, beta=_beta, dim=dim, spin=spin)
    println(UEG.short(para))
    kF = para.kF

    K = MCIntegration.Continuous(0.0, 200.0 * kF)
    result = integrate(integrand; userdata=para, var=K, neval=neval, parallel=:nothread)

    println(result)
    if isnothing(result) == false
        avg, std = result.mean[1] * spin, result.stdev[1] * spin
        open(filename, "a+") do io
            writedlm(io, [[_rs, _beta, spinPolarPara, avg, std]])
        end
    end
end