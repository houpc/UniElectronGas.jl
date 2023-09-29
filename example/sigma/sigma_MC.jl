using ElectronLiquid
using CompositeGrids

dim = 2
rs = [1.0,]
# mass2 = [1.0, 2.0, 3.0]
mass2 = [1.0]
Fs = [-0.0,]
beta = [25.0]
order = [3,]
neval = 2e7
isDynamic = false
isFock = false
diagGenerate = :GV
# diagGenerate = :Parquet

spin = 1
spinPolarPara = 2 / spin - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]

# mission = :Z
# mission = :K
mission = ARGS[1]
println("mission: ", mission)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock, spin=spin)
    println(UEG.short(para))
    kF = para.kF

    if mission == "Z"
        ######### calcualte Z factor ######################
        kgrid = [kF,]
        ngrid = [-1, 0]
    elseif mission == "K"
        ######### calculate K dependence #####################
        Nk, korder = 4, 4
        minK = 0.2kF
        # kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.2kF], [kF,], Nk, minK, korder).grid
        # kgrid = kF .+ [-0.1, -0.05, -0.03, -0.01, -0.005, -0.001, 0, 0.001, 0.005, 0.01, 0.03, 0.05, 0.1] * kF
        kgrid = kF .+ [-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1] * kF
        ngrid = [0,]
    else
        error("unknown mission")
    end

    filename = "data_$(mission)_test.jld2"
    # filename = "data$(dim)d_$(mission).jld2"

    sigma, result = Sigma.MC(para; kgrid=kgrid, ngrid=ngrid, spinPolarPara=spinPolarPara,
        neval=neval, filename=filename,
        diagtype=diagGenerate)
end
