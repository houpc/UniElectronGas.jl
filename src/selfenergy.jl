function sigma(para; neval=1e6, ngrid=[-1, 0, 1], kgrid=[para.kF], filename=nothing, diagtype=:GV, isLayered2D=false)
    sigma, result = Sigma.MC(para; kgrid=kgrid, ngrid=ngrid, neval=neval, filename=filename, diagtype=diagtype, isLayered2D=isLayered2D)
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
        # for p in sort([k for k in keys(data)])
        rdata[p] = real(sigma[p][:, :])
        idata[p] = imag(sigma[p][:, :])
        println(p, " ", rdata[p])
    end
    return ngrid, kgrid, rdata, idata
end

function mu(data, nidx=1, kidx=1) #assume the data are calculated with ngrid = [0], kgrid = [para.kF]
    return real(data[nidx, kidx])
end

# function zfactor_inverse(data, β, nids::Union{Tuple{Int,Int},Vector{Int}}=(1, 2), kidx=1) #assume the data are calculated with ngrid = [-1, 0], kgrid = [para.kF]
#     return @. (imag(data[nids[2], kidx]) - imag(data[nids[1], kidx])) / (2π / β)
# end

function zfactor_inverse(data, β, nidx::Int=1, kidx=1) #assume the data are calculated with n = 0, kgrid = [para.kF]
    return @. imag(data[nidx, kidx]) / (π / β)
end


# function zfactor_inverse(data, data1, β, nids=[1, 2], kidx=1) #assume the data are calculated with ngrid = [-1, 0], kgrid = [para.kF]
#     # return @. (imag(data[nids[2], kidx]) - imag(data[nids[1], kidx])) / (2π / β)
#     return @. (imag(data1[1, kidx]) - imag(data[2, kidx])) / (2π / β)
# end

function meff_inverse(data, para, kgrid, i_dk=1, nidx=1)
    kF_label = searchsortedfirst(kgrid, para.kF)
    println("kF_label: $kF_label")
    # return @. (data[nidx, kF_label+i_dk] - data[nidx, kF_label-i_dk]) / (kgrid[kF_label+i_dk] - kgrid[kF_label-i_dk]) * para.me / para.kF
    return @. (data[nidx, kF_label+i_dk] - data[nidx, kF_label]) / (kgrid[kF_label+i_dk] - kgrid[kF_label]) * para.me / para.kF
end

function get_dzmu(para, datatuple; verbose=0,
    # function get_dzmu(para, datatuple, datatuple1; verbose=0,
    parafile="para_wn_1minus0.csv", root_dir=@__DIR__, isSave=false)
    ngrid, kgrid, data = datatuple
    # ngrid, kgrid, data1 = datatuple1
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
        # _zinv[p] = zfactor_inverse(val, para.β, [2, 3])
        # _zinv[p] = zfactor_inverse(val, data1[p], para.β, [2, 3])
    end

    println(para.order)
    println(_mu)
    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)

    println("dmu: ", dmu)
    println("dz: ", dz)
    println("dzinv: ", dzinv)

    ########### save to csv  #################
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
    # para1 = ParaMC(rs=para.rs, beta=40.0, Fs=para.Fs, order=5, mass2=para.mass2, isDynamic=para.isDynamic, dim=para.dim, spin=para.spin)
    # _mu, _zinv = CounterTerm.getSigma(para1, parafile=parafile, root_dir=root_dir)
    _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)

    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)
    println("dmu: ", dmu)

    for (p, val) in rdata
        p == (6, 0, 0) && println("(6,0,0): ", val)
        p == (5, 0, 1) && println("(5,0,1): ", val)
        p == (5, 1, 0) && println("(5,1,0): ", val * dmu[1])
        p == (4, 1, 1) && println("(4,1,1): ", val * dmu[1])
        p == (4, 2, 0) && println("(4,2,0): ", val * dmu[1]^2)
        p == (3, 3, 0) && println("(3,3,0): ", val * dmu[1]^3)
        p == (2, 4, 0) && println("(2,4,0): ", val * dmu[1]^4)
        p == (1, 5, 0) && println("(1,5,0): ", val * dmu[1]^5)
        p == (4, 1, 0) && println("(4,1,0): ", val * dmu[2])
        p == (3, 2, 0) && println("(3,2,0): ", val * 2 * dmu[1] * dmu[2])
        p == (2, 3, 0) && println("(2,3,0): ", val * 3 * dmu[1]^2 * dmu[2])
        p == (2, 2, 0) && println("(2,2,0): ", val * (2 * dmu[1] * dmu[3] + dmu[2]^2))
        p == (3, 1, 0) && println("(3,1,0): ", val * dmu[3])
        p == (1, 1, 0) && println("(1,1,0): ", val * dmu[5])
    end
    println()
    for (p, val) in rdata
        p == (5, 0, 0) && println("(5,0,0): ", val)
        p == (4, 0, 1) && println("(4,0,1): ", val)
        p == (4, 1, 0) && println("(4,1,0): ", val * dmu[1])
        p == (3, 1, 1) && println("(3,1,1): ", val * dmu[1])
        p == (3, 2, 0) && println("(3,2,0): ", val * dmu[1]^2)
        p == (2, 3, 0) && println("(2,3,0): ", val * dmu[1]^3)
        p == (1, 4, 0) && println("(1,4,0): ", val * dmu[1]^4)
        p == (3, 1, 0) && println("(3,1,0): ", val * dmu[2])
        p == (2, 2, 0) && println("(2,2,0): ", val * 2 * dmu[1] * dmu[2])
        p == (1, 3, 0) && println("(1,3,0): ", val * 3 * dmu[1]^2 * dmu[2])
        p == (1, 2, 0) && println("(1,2,0): ", val * (2 * dmu[1] * dmu[3] + dmu[2]^2))
        p == (2, 1, 0) && println("(2,1,0): ", val * dmu[3])
        p == (1, 1, 0) && println("(1,1,0): ", val * dmu[4])
    end
    println()
    for (p, val) in rdata
        p == (4, 0, 0) && println("(4,0,0): ", val)
        p == (3, 0, 1) && println("(3,0,1): ", val)
        p == (3, 1, 0) && println("(3,1,0): ", val * dmu[1])
        p == (2, 1, 1) && println("(2,1,1): ", val * dmu[1])
        p == (2, 2, 0) && println("(2,2,0): ", val * dmu[1]^2)
        p == (1, 3, 0) && println("(1,3,0): ", val * dmu[1]^3)
        p == (2, 1, 0) && println("(2,1,0): ", val * dmu[2])
        p == (1, 2, 0) && println("(1,2,0): ", val * 2 * dmu[1] * dmu[2])
        p == (1, 1, 0) && println("(1,1,0): ", val * dmu[3])
    end

    rSw_k = CounterTerm.chemicalpotential_renormalization(para.order, rdata, dmu)
    # rSw_k = CounterTerm.chemicalpotential_renormalization(para.order, rdata, dmu, verbose=1)
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

function getMeff(para, rSigma, kgrid::Vector{Float64}, idx_dk::Int=1; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    order = para.order
    # ngrid, kgrid, rdata, idata = loaddata(para, filename)
    # _Meffinv = Dict()
    # for (p, val) in rdata
    #     _Meffinv[p] = -meff_inverse(val, para, kgrid, idx_dk)
    # end
    dMeffinv = []
    for rSig_eo in rSigma
        push!(dMeffinv, -meff_inverse(rSig_eo, para, kgrid, idx_dk))
    end
    println(dMeffinv)

    _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)

    dzinv, dmu, dz = CounterTerm.sigmaCT(order, _mu, _zinv)
    println("dmu: ", dmu)

    # dMeffinv, dmu, _dMeff = CounterTerm.sigmaCT(order, _mu, _Meffinv)
    # println("dMeffinv: ", dMeffinv)

    # println("_dMeff: ", _dMeff)
    zinv = Taylor1([1.0, dzinv...], order)

    # Meff = zinv * Taylor1([1.0, _dMeff...], order)
    _Meffinv = Taylor1([1.0, dMeffinv...], order)
    println("dReΣ / dk = ", _Meffinv)
    println(1.0 / _Meffinv)
    dMeff = [TaylorSeries.getcoeff(1.0 / _Meffinv, o) for o in 0:order]
    println(accumulate(+, dMeff))

    Meff = zinv / _Meffinv
    dMeff = [TaylorSeries.getcoeff(Meff, o) for o in 1:order]

    sumMeff = accumulate(+, dMeff)
    return @. 1.0 + sumMeff
end

function getMeff(para, rSigma, kgrid::Vector{Float64}; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    order, kF = para.order, para.kF

    # idxs = [2, 3, 4, 5]
    # fit ReΣ(k,ω0) = a*(k-kF) + b*(k-kF)^2
    dMeffinv = []
    x = kgrid
    # x = kgrid[idxs]
    @. model(x, p) = p[1] + p[2] * (x - kF) + p[3] * (x - kF)^2
    println("ReΣ :", rSigma)
    fit_parameters = []
    for o in 1:order
        ydat = rSigma[o]
        # ydat = rSigma[o][idxs]
        y = vec(Measurements.value.(ydat))
        wt = vec(1 ./ Measurements.uncertainty.(ydat) .^ 2)
        fit = curve_fit(model, x, y, wt, [0.1, -0.1, 0.0])
        push!(dMeffinv, -Measurements.measurement(coef(fit)[2], stderror(fit)[2]) * para.me / kF)
        push!(fit_parameters, [coef(fit), stderror(fit)])
    end

    # para1 = ParaMC(rs=para.rs, beta=40.0, Fs=para.Fs, order=para.order, mass2=para.mass2, isDynamic=para.isDynamic, dim=para.dim, spin=para.spin)
    # _mu, _zinv = CounterTerm.getSigma(para1, parafile=parafile, root_dir=root_dir)
    _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)
    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)
    zinv = Taylor1([1.0, dzinv...], order)
    println("1 / z = ", zinv)

    _Meffinv = Taylor1([1.0, dMeffinv...], order)
    println("dReΣ / dk = ", _Meffinv)
    Meff = zinv / _Meffinv
    dMeff = [TaylorSeries.getcoeff(Meff, o) for o in 1:order]
    println("δm* / m = ", dMeff)

    sumMeff = accumulate(+, dMeff)
    return @. 1.0 + sumMeff, fit_parameters
    # return @. 1.0 + sumMeff
end

function getMeff(para, rSigmadk; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    order, kF = para.order, para.kF

    dMeffinv = -rSigmadk .* para.me / kF
    # para1 = ParaMC(rs=para.rs, beta=40.0, Fs=para.Fs, order=para.order, mass2=para.mass2, isDynamic=para.isDynamic, dim=para.dim, spin=para.spin)
    # _mu, _zinv = CounterTerm.getSigma(para1, parafile=parafile, root_dir=root_dir)
    _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)
    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)
    zinv = Taylor1([1.0, dzinv...], order)
    println("1 / z = ", zinv)

    _Meffinv = Taylor1([1.0, dMeffinv...], order)
    println("dReΣ / dk = ", _Meffinv)
    Meff = zinv / _Meffinv
    dMeff = [TaylorSeries.getcoeff(Meff, o) for o in 1:order]
    println("δm* / m = ", dMeff)

    sumMeff = accumulate(+, dMeff)
    return @. 1.0 + sumMeff
    # return @. 1.0 + sumMeff
end