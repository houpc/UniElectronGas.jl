using ElectronGas
using ElectronGas.CompositeGrids
using ElectronGas.Interaction
using ElectronGas.Parameters
const EPS = 1e-16


const rs = [1.0, 2.0, 3.0, 4.0, 5.0]
const Fs = -[0.223, 0.380, 0.516, 0.639, 0.752]
const n = 0

function dielectric(q, n, param; Fs=0.0)
    @unpack spin = param

    if abs(q) < EPS
        q = EPS
    end

    Π::Float64 = spin * Polarization.Polarization0_3dZeroTemp(q, n, param)

    fs = Fs
    Vinvs, Vinva = Interaction.coulombinv(q, param)
    Ks = Vinvs / (Vinvs - (1 + fs * Vinvs) * Π)
    return Ks
end

function ver3AA(rs, Fs)
    param = Parameter.rydbergUnit(0.001, rs)
    xgrid = ElectronGas.CompositeGrids.CompositeG.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 32, 1e-16, 8)
    qs = [sqrt(param.kF^2 * 2 - 2 * x * param.kF^2) for x in xgrid]
    data = [dielectric(q, n, param; Fs=Fs) for q in qs]
    return ElectronGas.CompositeGrids.Interp.integrate1D(data, xgrid)
end

println([ver3AA(rs[i], Fs[i]) for i in 1:length(rs)])