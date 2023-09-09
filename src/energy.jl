function getEnergy(para, filename; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    # if para.order < 4
    #     para1 = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=para.Fs, Fa=-0.0, order=4, dim=para.dim,
    #         mass2=para.mass2, isDynamic=para.isDynamic, isFock=para.isFock)
    # end
    para1 = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=para.Fs, Fa=-0.0, order=3, dim=para.dim,
        mass2=para.mass2, isDynamic=para.isDynamic, isFock=para.isFock)
    _mu, _zinv = CounterTerm.getSigma(para1, parafile=parafile, root_dir=root_dir)
    dzinv, dmu, dz = CounterTerm.sigmaCT(para1.order, _mu, _zinv)

    key = UEG.short(para)
    f = jldopen(filename, "r")
    data = f[key][1]

    println(data)
    return CounterTerm.chemicalpotential_renormalization(para.order, data, dmu)
end