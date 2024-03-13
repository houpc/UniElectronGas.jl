using ElectronGas, ElectronLiquid, Parameters
using GreenFunc, CompositeGrids
using PyCall

@pyimport numpy as np   # for saving/loading numpy data

@inline function get_Fs(rs)
    return get_Fs_PW(rs)
    # if rs < 1.0 || rs > 5.0
    #     return get_Fs_DMC(rs)
    # else
    #     return get_Fs_PW(rs)
    # end
end

"""
Get the symmetric l=0 Fermi-liquid parameter F⁰ₛ from DMC data of 
Moroni, Ceperley & Senatore (1995) [Phys. Rev. Lett. 75, 689].
"""
@inline function get_Fs_DMC(rs)
    return error("Not yet implemented!")
end

"""
Get the symmetric l=0 Fermi-liquid parameter F⁰ₛ via interpolation of the 
compressibility ratio data of Perdew & Wang (1992) [Phys. Rev. B 45, 13244].
"""
@inline function get_Fs_PW(rs)
    # if rs < 1.0 || rs > 5.0
    #     @warn "The Perdew-Wang interpolation for Fs may " *
    #           "be inaccurate outside the metallic regime!"
    # end
    kappa0_over_kappa = 1.0025 - 0.1721rs - 0.0036rs^2
    # F⁰ₛ = κ₀/κ - 1
    return kappa0_over_kappa - 1.0
end

"""
Get the antisymmetric l=0 Fermi-liquid parameter F⁰ₐ via interpolation of the 
susceptibility ratio data (c.f. Kukkonen & Chen, 2021)
"""
@inline function get_Fa(rs)
    chi0_over_chi = 0.9821 - 0.1232rs + 0.0091rs^2
    # F⁰ₐ = χ₀/χ - 1
    return chi0_over_chi - 1.0
end

# exchange interaction (Ws + Wa \sigma\sigma)_ex to a direct interaction Ws'+Wa' \sigma\sigma 
# # exchange S/A interaction projected to the spin-symmetric and antisymmetric parts
function exchange2direct(Wse, Wae)
    Ws = (Wse + 3 * Wae) / 2
    Wa = (Wse - Wae) / 2
    return Ws, Wa
end

function main()
    # System parameters
    dim = 2
    # dim = 3
    mass2 = 1e-8
    beta = 1000.0

    # Output directory
    mkpath("results/finalized_meff_results/$(dim)d")
    dir = joinpath(@__DIR__, "results/finalized_meff_results/$(dim)d")

    rslist = collect(LinRange(0.01, 1.1, 41))  # 2D
    # rslist = collect(LinRange(0.002, 1.1, 111))  # 2D
    # rslist = collect(LinRange(0.01, 10.0, 101))  # 3D

    # NOTE: int_type ∈ [:ko_const, :ko_takada_plus, :ko_takada, :ko_moroni, :ko_simion_giuliani] 
    # NOTE: KO interaction using G+ is currently only available in 3D
    int_type = :ko_const
    # int_type = :ko_moroni
    # int_type = :ko_takada_plus
    @assert int_type ∈ [:ko_const, :ko_takada_plus, :ko_moroni]

    if int_type == :ko_const
        ko_dirstr = "const"
    elseif int_type == :ko_takada_plus
        ko_dirstr = "takada"
    elseif int_type == :ko_moroni
        ko_dirstr = "simion_giuliani"
    end

    F1_rpa_vs_rs = []
    F1_rpa_fl_vs_rs = []
    for rs in rslist
        param = Parameter.rydbergUnit(1.0 / beta, rs, dim)
        kF = param.kF

        # Small params
        Euv, rtol = 100 * param.EF, 1e-11
        # Euv, rtol = 1000 * param.EF, 1e-11
        # Euv, rtol = 1000 * param.EF, 1e-14

        # # DLR cutoffs
        # Euv, rtol = 1000 * param.EF, 1e-11
        # # maxK, minK = 20kF, 1e-8kF
        # # Nk, order = 12, 8

        if dim == 3
            nu_grid = CompositeGrid.LogDensedGrid(
                :gauss,
                [-1.0, 1.0],
                [-1.0, 0.0, 1.0],
                20,
                1e-6,
                12,
            )
            # nu_grid = CompositeGrid.LogDensedGrid(
            #     :gauss,
            #     [-1.0, 1.0],
            #     [-1.0, 0.0, 1.0],
            #     16,
            #     1e-6,
            #     32,
            # )
            # nu_grid = cos.(CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 32))
            # 3D: |k - k'| = kF sqrt(2(1 - ν)), ν = cosθ
            costhetas = nu_grid.grid
            # costhetas = nu_grid
        else
            # dt = 1e-7
            # dt = 2.177349558186689e-7
            # theta_grid = CompositeGrid.LogDensedGrid(
            #     :gauss,
            #     [dt, 2π - dt],
            #     [dt, π, 2π - dt],
            #     20,
            #     1e-6,
            #     12,
            # )
            theta_grid = CompositeGrid.LogDensedGrid(:gauss, [0.0, 2π], [0.0, 2π], 16, 0.001, 32)
            # 2D: |k - k'| = kF sqrt(2(1 - cosθ))
            costhetas = cos.(theta_grid.grid)
            sinthetas = sin.(theta_grid.grid)
        end

        # k_m_kps = @. kF * sqrt(2 * (1 - costhetas))
        k_m_kps = @. sqrt(2 * kF^2 - 2 * costhetas * kF^2)

        if dim == 3
            # Get Fermi liquid parameter F⁰ₛ(rs) from Perdew-Wang fit
            rs = round(param.rs; sigdigits=13)
            Fs = get_Fs(rs)
            Fa = get_Fa(rs)
            if int_type in [:ko_const_p, :ko_const_pm]
                _int_type = :ko_const
                if int_type == :ko_const_pm
                    println("Fermi liquid parameters at rs = $(rs): Fs = $Fs, Fa = $Fa")
                end
            else
                _int_type = int_type
            end

            # Either use constant Fs from P&W, q-dependent Takada ansatz, or Corradini fit to Moroni DMC data
            if int_type == :ko_const
                landaufunc = Interaction.landauParameterConst
            elseif int_type == :ko_takada_plus
                landaufunc = Interaction.landauParameterTakadaPlus
            elseif int_type == :ko_moroni
                landaufunc = Interaction.landauParameterMoroni
            end
        end

        # Get W0_{+}(p, iωₙ=0) * Vinv (spin-symmetric static exchange part of RPA interaction divided by V)
        W0_dyn, W0_inst = Interaction.RPAwrapped(
            Euv,
            rtol,
            k_m_kps,
            param;
            regular=false,
        )
        # W0_inst = 1 / V, W0_dyn = Wtilde, so W = (1 / W0_inst + W0_dyn) = V + Wtilde
        W0_plus_static = @. (1 / W0_inst[1, 1, :] + W0_dyn[1, 1, :])
        @assert maximum(imag(W0_plus_static)) ≤ 1e-10
        # println(maximum(imag(W0_plus_static)))

        if dim == 3
            # Get R_{+}(p, iωₙ=0) * Vinv (spin-symmetric static exchange part of KO interaction divided by V)
            R_dyn, R_inst = Interaction.KOwrapped(
                Euv,
                rtol,
                k_m_kps,
                param;
                regular=false,
                int_type=_int_type,
                landaufunc=landaufunc,
                Fs=_int_type == :ko_const ? Fs : -0.0,   # NOTE: NEFT uses opposite sign convention for F!
                Fa=int_type == :ko_const_pm ? Fa : -0.0,  # NOTE: NEFT uses opposite sign convention for F!
                # Fs=_int_type == :ko_const ? -Fs : -0.0,   # NOTE: NEFT uses opposite sign convention for F!
                # Fa=int_type == :ko_const_pm ? -Fa : -0.0,  # NOTE: NEFT uses opposite sign convention for F!
            )
            # R_inst = 1 / V, R_dyn = Rtilde, so R = (1 / R_inst + R_dyn) = V + Rtilde
            R_plus_static = @. (1 / R_inst[1, 1, :] + R_dyn[1, 1, :])
            @assert maximum(imag(R_plus_static)) ≤ 1e-10
        end

        # The RPA and RPA+FL integrands are ν W0_{+}(p, iωₙ=0) and ν R_{+}(p, iωₙ=0), respectively
        if dim == 3
            F1p_rpa_integrand = @. -param.NF * costhetas * real(W0_plus_static)
            F1p_rpa_fl_integrand = @. -param.NF * costhetas * real(R_plus_static)
        else
            F1p_rpa_integrand = @. -param.NF * costhetas * real(W0_plus_static)
        end

        # Perform angular integrations ν ∈ [-1, 1]
        if dim == 3
            F1p_rpa = 0.25 * CompositeGrids.Interp.integrate1D(F1p_rpa_integrand, nu_grid)
            F1p_rpa_fl = 0.25 * CompositeGrids.Interp.integrate1D(F1p_rpa_fl_integrand, nu_grid)
        else
            F1p_rpa = (1 / 4π) * CompositeGrids.Interp.integrate1D(F1p_rpa_integrand, theta_grid)
        end

        # Use ElectronLiquid.jl to compute the same quantities
        p_rpa = ParaMC(; rs=rs, beta=beta, dim=dim, Fs=-0.0, order=1, mass2=mass2)
        if dim == 3
            p_rpa_fl = ParaMC(; rs=rs, beta=beta, dim=dim, Fs=-Fs, order=1, mass2=mass2)  # NOTE: ElectronLiquid uses opposite convention for Fs compared to ElectronGas!
        end

        # NOTE: NEFT uses opposite sign convention for F!
        if dim == 3
            F1p_rpa_ueg, _ = -0.5 .* Ver4.projected_exchange_interaction(1, p_rpa, Ver4.exchange_interaction)
            F1p_rpa_fl_ueg, _ = -0.5 .* Ver4.projected_exchange_interaction(1, p_rpa_fl, Ver4.exchange_interaction)
        else
            F1p_rpa_ueg, _ = -(1 / 2π) .* Ver4.projected_exchange_interaction(1, p_rpa, Ver4.exchange_interaction)
        end

        if dim == 3
            println("\nrs = $rs:" * "\nF^{(RPA)+}_1 = $F1p_rpa" * "\nF^{(RPA+FL)+}_1 = $F1p_rpa_fl")
            println("\nElectronLiquid:" * "\nF^{(RPA)+}_1 = $F1p_rpa_ueg" * "\nF^{(RPA+FL)+}_1 = $F1p_rpa_fl_ueg")
        else
            println("\nrs = $rs:" * "\nF^{(RPA)+}_1 = $F1p_rpa")
            println("\nElectronLiquid:" * "\nF^{(RPA)+}_1 = $F1p_rpa_ueg")
        end

        # push!(F1_rpa_vs_rs, F1p_rpa)
        push!(F1_rpa_vs_rs, F1p_rpa_ueg)
        if dim == 3
            # push!(F1_rpa_fl_vs_rs, F1p_rpa_fl)
            push!(F1_rpa_fl_vs_rs, F1p_rpa_fl_ueg)
        end
    end

    # Add points at rs = 0
    pushfirst!(rslist, 0.0)
    pushfirst!(F1_rpa_vs_rs, 0.0)
    if dim == 3
        pushfirst!(F1_rpa_fl_vs_rs, 0.0)
    end

    # Save data
    np.savez(joinpath(dir, "rpa/meff_$(dim)d_tree_level_G0W0_test.npz"), rslist=rslist, mefflist=(1 .+ F1_rpa_vs_rs))
    if dim == 3
        np.savez(joinpath(dir, "$(ko_dirstr)/meff_$(dim)d_tree_level_G0Wp_test.npz"), rslist=rslist, mefflist=(1 .+ F1_rpa_fl_vs_rs))
    end
end

main()
