using ElectronLiquid, FeynmanDiagram
using CompositeGrids

dim = 3
rs = [1.0,]
mass2 = [1.0]
Fs = [-0.0,]
beta = [25.0]
order = [3,]
neval = 2e7
isDynamic = false
isFock = false
diagGenerate = :GV
# diagGenerate = :Parquet
response = ChargeCharge
# response = SpinSpin

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    println(UEG.short(para))
    # para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock, spin=1)
    kF = para.kF

    ######### calculate (K,n) dependence #####################
    Nk, korder = 4, 4
    minK = 0.2kF
    # kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.2kF], [kF,], Nk, minK, korder).grid
    kgrid = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] * kF
    ngrid = [0,]

    filename = "data$(dim)d_polar_test.jld2"

    sigma, result = Polarization.MC(para; kgrid=kgrid, ngrid=ngrid, neval=neval, filename=filename,
        response=response, diagtype=diagGenerate)

end
