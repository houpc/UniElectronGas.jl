using ElectronLiquid, FeynmanDiagram, JLD2

function ver3_static(para::ParaMC, diagram, anglegrid;
    filename=nothing,
    kamp=para.kF,
    kamp2=kamp,
    kwargs...)
    kin = [kamp, 0, 0]
    kouts = [[kamp2 * cos(θ), kamp2 * sin(θ), 0] for θ in anglegrid]
    qout = [(kin .- k) for k in kouts]

    ver3, result = Ver3.KW(para, diagram;
        kin=[kin,],
        qout=qout,
        kwargs...)

    if (isnothing(result) == false)
        if (isnothing(filename) == false)
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                println(key)
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kamp, kamp2, anglegrid, ver3)
            end
        end
        return ver3, result
    else
        return nothing, nothing
    end

end

function ver3_angle(para::ParaMC, diagram;
    filename=nothing,
    kamp=para.kF,
    kamp2=kamp,
    kwargs...)
    ver3, result = Ver3.AA(para, diagram;
        kamp=[kamp,],
        kamp2=[kamp2,],
        kwargs...)
    if (isnothing(result) == false)
        if (isnothing(filename) == false)
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                println(key)
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kamp, kamp2, [0,], ver3)
            end
        end
        return ver3, result
    else
        return nothing, nothing
    end

end

function ver3_KW(para::ParaMC, diagram;
    kin=[getK(para.kF, para.dim, 1),],
    nkin=[0,],
    qout=[getK(0.0, para.dim, 1),],
    nqout=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    filename=nothing,
    kwargs...)
    ver3, result = Ver3.KW(para, diagram;
        kin=kin,
        nkin=nkin,
        qout=qout, nqout=nqout,
        neval=neval,
        print=print,
        alpha=alpha,
        config=config,
        kwargs...)

    if (isnothing(result) == false)
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kin, qout, nqout, ver3)
                return ver3, result
            end
        end
    else
        return nothing, nothing
    end

end
