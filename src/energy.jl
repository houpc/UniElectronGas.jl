function getEnergy(para, filename, E0=nothing; parafile="para_wn_1minus0.csv", root_dir=@__DIR__)
    # if para.order < 4
    #     para1 = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=para.Fs, Fa=-0.0, order=4, dim=para.dim,
    #         mass2=para.mass2, isDynamic=para.isDynamic, isFock=para.isFock)
    # else
    #     para1 = para
    # end
    _mu, _zinv = CounterTerm.getSigma(para, parafile=parafile, root_dir=root_dir)
    dzinv, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _zinv)

    density, order = para.basic.n, para.order
    key = UEG.short(para)
    f = jldopen(filename, "r")
    data = f[key][1]

    println("dmu: ", dmu)
    dF_eachorder = CounterTerm.chemicalpotential_renormalization(order, data, dmu)
    println("dF: ", dF_eachorder)

    # _dF0 = order0_renormalize(order, data, dmu)
    # dF0 = Taylor1([density, _dF0...])
    # dF0 = dF0 * Taylor1(dmu)
    # dF0_eachorder = [getcoeff(dF0, o) for o in 0:order-1]

    dF0_eachorder = order0_renormalization(para.order, data, dmu, density)
    println("dF0: ", dF0_eachorder)

    if isnothing(E0)
        E_eachorder = [data[(0, 0, 0)],]
    else
        E_eachorder = [measurement(E0[1], E0[2]),]
    end
    # append!(E_eachorder, dF_eachorder .+ dF0_eachorder)
    append!(E_eachorder, dF_eachorder)

    println("free energy/V from each order: ", E_eachorder)
    F = accumulate(+, E_eachorder) ./ density
    println("free energy/N at each order (Ha): ", F ./ 2)
    return F
end


function order0_renormalization(order, data, δμ, density=1.0)
    @assert order <= 5 "Order $order hasn't been implemented!"
    @assert length(δμ) + 1 >= order
    d = CounterTerm.mergeInteraction(data)
    sample = collect(values(d))[1]
    z = [zero(sample) for i in 1:order]

    if order >= 1
        z[1] = density * δμ[1]
    end
    if order >= 2
        z[2] = density * δμ[2] + d[(0, 1)] * δμ[1]^2
    end
    if order >= 3
        z[3] = density * δμ[3] + d[(0, 1)] * 2 * δμ[1] * δμ[2] + d[(0, 2)] * δμ[1]^3
    end
    if order >= 4
        z[4] = density * δμ[4] +
               d[(0, 1)] * (δμ[2]^2 + 2 * δμ[1] * δμ[3]) +
               d[(0, 2)] * 3 * δμ[1]^2 * δμ[2] +
               d[(0, 3)] * δμ[1]^4
    end
    if order >= 5
        z[5] = density * δμ[5] +
               d[(0, 1)] * (2 * δμ[1] * δμ[4] + 2 * δμ[2] * δμ[3]) +
               d[(0, 2)] * (3 * δμ[1]^2 * δμ[3] + 3 * δμ[1] * δμ[2]^2) +
               d[(0, 3)] * 4 * δμ[1]^3 * δμ[2] +
               d[(0, 4)] * δμ[1]^5
    end

    return z
end

function order0_renormalize(order, data, δμ)
    @assert order <= 5 "Order $order hasn't been implemented!"
    @assert length(δμ) + 1 >= order
    d = CounterTerm.mergeInteraction(data)
    sample = collect(values(d))[1]
    z = [zero(sample) for i in 1:order]

    if order >= 1
        z[1] = d[(0, 1)] * δμ[1]
    end
    if order >= 2
        z[2] = d[(0, 1)] * δμ[2] + d[(0, 2)] * δμ[1]^2
    end
    if order >= 3
        z[3] = d[(0, 1)] * δμ[3] +
               d[(0, 2)] * 2 * δμ[1] * δμ[2] +
               d[(0, 3)] * δμ[1]^3
    end
    if order >= 4
        z[4] = d[(0, 1)] * δμ[4] +
               d[(0, 2)] * (2 * δμ[1] * δμ[3] + δμ[2]^2) +
               d[(0, 3)] * 3 * δμ[1]^2 * δμ[2] +
               d[(0, 4)] * δμ[1]^4
    end

    return z
end