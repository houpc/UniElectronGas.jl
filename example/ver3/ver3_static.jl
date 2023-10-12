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

    if isnothing(filename) == false
        jldopen(filename, "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (kamp, kamp2, anglegrid, ver3)
        end
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