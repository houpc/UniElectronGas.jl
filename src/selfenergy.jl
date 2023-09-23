function sigma(para; neval=1e6, ngrid=[-1, 0, 1], kgrid=[para.kF], filename=nothing, diagtype=:GV)
    sigma, result = Sigma.MC(para; kgrid=kgrid, ngrid=ngrid, neval=neval, filename=filename, diagtype=diagtype)
    return sigma, result
end

function loaddata(para, FileName)
    key = UEG.short(para)
    f = jldopen(FileName, "r")

    p = ParaMC(key)
    ngrid, kgrid, sigma = f[key]
    _partition = UEG.partition(para.order)
    rdata, idata = Dict(), Dict()
    for p in _partition
        rdata[p] = real(sigma[p][:, :])
        idata[p] = imag(sigma[p][:, :])
        println(p, " ", rdata[p])
    end
    return ngrid, kgrid, rdata, idata
end

function mu(data, nidx=2, kidx=1) #assume the data are calculated with ngrid = [-1, 0], kgrid = [para.kF]
    return real(data[nidx, kidx])
end

function zfactor_inverse(data, β, nids=[1, 2], kidx=1) #assume the data are calculated with ngrid = [-1, 0], kgrid = [para.kF]
    return @. (imag(data[nids[2], kidx]) - imag(data[nids[1], kidx])) / (2π / β)
end

function meff_inverse(data, para, kgrid, i_dk=1, nidx=1)
    kF_label = searchsortedfirst(kgrid, para.kF)
    return @. (data[nidx, kF_label+i_dk] - data[nidx, kF_label-i_dk]) / (kgrid[kF_label+i_dk] - kgrid[kF_label-i_dk]) * para.me / para.kF
end

function get_dzmu(para, datatuple; verbose=0,
    parafile="para_wn_1minus0.csv", root_dir=@__DIR__, isSave=false)
    ngrid, kgrid, data = datatuple
    printstyled(UEG.short(para), color=:yellow)
    println()

    if verbose > 0
        for p in sort([k for k in keys(data)])
            println("$p: μ = $(mu(data[p]))   z_inverse = $(zfactor_inverse(data[p], para.β))")
        end
    end

    _mu = Dict()
    for (p, val) in data
        _mu[p] = mu(val)
    end
    _zinv = Dict()
    for (p, val) in data
        _zinv[p] = zfactor_inverse(val, para.β)
    end

    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)

    println("dmu: ", dmu)
    println("dz: ", dz)
    println("dzinv: ", dzinv)

    ############# save to csv  #################
    if isSave
        df = CounterTerm.fromFile(parafile, root_dir=root_dir)
        for P in keys(data)
            paraid = UEG.paraid(para)
            df = CounterTerm.appendDict(df, paraid, Dict("partition" => P, "μ" => _mu[P].val, "μ.err" => _mu[P].err, "Σw" => _zinv[P].val, "Σw.err" => _zinv[P].err); replace=true)
        end
        CounterTerm.toFile(df, parafile, root_dir=root_dir)
    end

    return dz, dzinv, dmu
end

function getdmu(para; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)
    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)
    return dmu
end

function getSigma(para, filename=filename; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    ngrid, kgrid, rdata, idata = loaddata(para, filename)
    para1 = ParaMC(rs=para.rs, beta=40.0, Fs=para.Fs, order=5, mass2=para.mass2, isDynamic=para.isDynamic, dim=para.dim, spin=para.spin)
    _mu, _zinv = CounterTerm.getSigma(para1, parafile=parafile, root_dir=root_dir)
    # _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)

    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)
    rSw_k = CounterTerm.chemicalpotential_renormalization(para.order, rdata, dmu, verbose=1)
    iSw_k = CounterTerm.chemicalpotential_renormalization(para.order, idata, dmu)

    return ngrid, kgrid, rSw_k, iSw_k
end

function getZfactor(para; parafile="para_wn_1minus0.csv", root_dir=@__DIR__, isRenorm=true)
    _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)
    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)

    if isRenorm
        sumzinv = accumulate(+, dzinv)
        return @. 1.0 / (1.0 + sumzinv)
    else
        sumz = accumulate(+, dz)
        return @. 1.0 + sumz
    end
end

function getMeff(para, filename, idx_dk::Int=1; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    order = para.order
    ngrid, kgrid, rdata, idata = loaddata(para, filename)
    _Meffinv = Dict()
    for (p, val) in rdata
        _Meffinv[p] = -meff_inverse(val, para, kgrid, idx_dk)
    end

    _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)

    dzinv, dmu, dz = CounterTerm.sigmaCT(order, _mu, _zinv)
    println("dmu: ", dmu)

    dMeffinv, dmu, _dMeff = CounterTerm.sigmaCT(order, _mu, _Meffinv)
    println("dMeffinv: ", dMeffinv)

    println("_dMeff: ", _dMeff)
    zinv = Taylor1([1.0, dzinv...], order)

    Meff = zinv * Taylor1([1.0, _dMeff...], order)
    dMeff = [getcoeff(Meff, o) for o in 1:order]

    sumMeff = accumulate(+, dMeff)
    return @. 1.0 + sumMeff
end

function getMeff(para, rSigma, kgrid::Vector{Float64}; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    order, kF = para.order, para.kF

    # fit ReΣ(k,ω0) = a*(k-kF) + b*(k-kF)^2
    dMeffinv = []
    x = kgrid
    @. model(x, p) = p[1] + p[2] * (x - kF) + p[3] * (x - kF)^2
    println(rSigma)
    fit_parameters = []
    for o in 1:order
        y = vec(Measurements.value.(rSigma[o]))
        wt = vec(1 ./ Measurements.uncertainty.(rSigma[o]) .^ 2)
        fit = curve_fit(model, x, y, wt, [1.0, -0.1, 0.0])
        push!(dMeffinv, -Measurements.measurement(coef(fit)[2], stderror(fit)[2]) * para.me / kF)
        push!(fit_parameters, [coef(fit), stderror(fit)])
    end

    para1 = ParaMC(rs=para.rs, beta=40.0, Fs=para.Fs, order=para.order, mass2=para.mass2, isDynamic=para.isDynamic, dim=para.dim, spin=para.spin)
    _mu, _zinv = CounterTerm.getSigma(para1, parafile=parafile, root_dir=root_dir)
    # _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)
    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)
    zinv = Taylor1([1.0, dzinv...], order)
    println(zinv)

    _Meffinv = Taylor1([1.0, dMeffinv...], order)
    println(_Meffinv)
    Meff = zinv / _Meffinv
    dMeff = [getcoeff(Meff, o) for o in 1:order]

    sumMeff = accumulate(+, dMeff)
    return @. 1.0 + sumMeff, fit_parameters
    # return @. 1.0 + sumMeff
end