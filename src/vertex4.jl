function ver4_PH(para; kamp=[para.kF,], kamp2=kamp, n=[-1, 0, 0, -1], ell=[0,],
    neval=1e6, filename::Union{String,Nothing}=nothing,
    partition=UEG.partition(para.order)
)
    ver4, result = Ver4.MC_PH(para, kamp=kamp, kamp2=kamp2, n=n, l=ell, neval=neval,
        filename=filename, partition=partition, fileter=[NoHartree])
    return ver4, result
end

function getVer4PHl(para, filename; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    if para.order < 5
        para1 = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=para.Fs, Fa=-0.0, order=5, dim=para.dim,
            mass2=para.mass2, isDynamic=para.isDynamic, isFock=para.isFock)
        _mu, _zinv = CounterTerm.getSigma(para1, parafile=parafile, root_dir=root_dir)
        dzinv, dmu, dz = CounterTerm.sigmaCT(para1.order, _mu, _zinv)
    else
        _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)
        dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)
    end

    vuu, vud = ver4PHl_renormalization(para, filename, dz, dmu)

    vuu = accumulate(+, vuu)
    vud = accumulate(+, vud)
    return (vuu + vud) / 2.0, (vuu - vud) / 2.0
end

function ver4PHl_renormalization(para, filename, dz, dmu)
    # println("read Fs = $Fs from $filename")
    kF = para.kF
    f = jldopen(filename, "r")
    # z1 = zeros(Measurement{Float64}, length(Fs), length(Λgrid))

    vuu = Dict()
    vud = Dict()

    key = UEG.short(para)
    kgrid, n, l, ver4 = f[key]

    for p in keys(ver4)
        println(p)
        if haskey(vuu, p) == false
            vuu[p] = MeshArray(l, kgrid; dtype=Complex{Measurement{Float64}})
            vud[p] = MeshArray(l, kgrid; dtype=Complex{Measurement{Float64}})
        end
        vuu[p][:, :] = ver4[p][1, :, :]
        vud[p][:, :] = ver4[p][2, :, :]
    end

    vuu_renorm = [vuu[(1, 0, 0)],]
    # sample = collect(values(vuu))[1]
    # z = [zero(sample) for i in 1:order]
    append!(vuu_renorm, CounterTerm.chemicalpotential_renormalization(para.order - 1, vuu, dmu, offset=1))
    vuu_renorm = CounterTerm.z_renormalization(para.order, vuu_renorm, dz, 2) #left leg renormalization

    vud_renorm = [vud[(1, 0, 0)],]
    append!(vud_renorm, CounterTerm.chemicalpotential_renormalization(para.order - 1, vud, dmu, offset=1))
    vud_renorm = CounterTerm.z_renormalization(para.order, vud_renorm, dz, 2) #left leg renormalization

    # vuu = [vuu[(o, 0)] for o in 1:para.order]
    # vud = [vud[(o, 0)] for o in 1:para.order]
    return vuu_renorm, vud_renorm
end
