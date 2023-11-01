dim = 3 # dimension of the problem
rs = [5.0]
mass2 = [0.875]  # screening parameter
Fs = [-0.0]    # Fermi liquid parameter with zero angular momentum
beta = [40.0]   # inverse temperature beta = β*E_F 
order = [5]    # order of diagrams
neval = 1e11    # number of Monte Carlo samples
isDynamic = false # whether to use effective field theory with dynamic screening or not 
isFock = false # whether to use Fock renormalization or not

diagGenerate = :Parquet # :GV or :Parquet, algorithm to generate diagrams
# diagGenerate = :GV # :GV or :Parquet, algorithm to generate diagrams
isLayered2D = false # whether to use layered 2D system or not

spin = 2    # 2 for unpolarized, 1 for polarized
spinPolarPara = 2 / spin - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
ispolarized = spinPolarPara != 0.0

# Whether or not we are benchmarking the diagrams
benchmark = :Parquet
# benchmark = :GV_new
# benchmark = :GV_old
# benchmark = nothing
if benchmark == :Parquet
    data_directory = joinpath(@__DIR__, "sigma/data$(dim)d/parquet")
    res_directory = joinpath(@__DIR__, "sigma/benchmark/parquet")
    para_directory = res_directory
elseif benchmark == :GV_new
    data_directory = joinpath(@__DIR__, "sigma/data$(dim)d/gv_new")
    res_directory = joinpath(@__DIR__, "sigma/benchmark/gv_new")
    para_directory = res_directory
elseif benchmark == :GV_old
    data_directory = joinpath(@__DIR__, "sigma/data$(dim)d/gv_old")
    res_directory = joinpath(@__DIR__, "sigma/benchmark/gv_old")
    para_directory = res_directory
else
    data_directory = joinpath(@__DIR__, "sigma/data$(dim)d")
    res_directory = joinpath(@__DIR__, "sigma")
    para_directory = ""
end

if isLayered2D
    const parafilename = joinpath(para_directory, "para_wn_1minus0_layered2d.csv")
    const simga_z_filename = joinpath(data_directory, "data$(dim)d_Z_$(diagGenerate)_layered.jld2")
    const sigma_k_filename = joinpath(data_directory, "data$(dim)d_K_$(diagGenerate)_layered.jld2")
    if spin == 2
        const meff_filename = joinpath(res_directory, "meff_$(dim)d_$(diagGenerate)_layered.dat")
        const zfactor_filename = joinpath(res_directory, "zfactor_layered2d.dat")
    else
        const meff_filename = joinpath(res_directory, "meff_$(dim)d_$(diagGenerate)_layered_spin$(spin).dat")
        const zfactor_filename = joinpath(res_directory, "zfactor_layered2d_spin$(spin).dat")
    end
elseif ispolarized
    const parafilename = joinpath(para_directory, "para_wn_1minus0_GV_spin_polarized.csv")
    if spin == 2
        const meff_filename = joinpath(res_directory, "meff_$(dim)d_GV_spin_polarized.dat")
        const zfactor_filename = joinpath(res_directory, "zfactor_$(dim)d_GV_spin_polarized.dat")
    else
        const meff_filename = joinpath(res_directory, "meff_$(dim)d_spin$(spin)_GV_spin_polarized.dat")
        const zfactor_filename = joinpath(res_directory, "zfactor_$(dim)d_spin$(spin)_GV_spin_polarized.dat")
    end
    const sigma_z_filename = joinpath(data_directory, "data$(dim)d_Z_GV_spin_polarized.jld2")
    const sigma_k_filename = joinpath(data_directory, "data$(dim)d_K_GV_spin_polarized.jld2")
else
    const parafilename = joinpath(para_directory, "para_wn_1minus0.csv")
    if spin == 2
        const meff_filename = joinpath(res_directory, "meff_$(dim)d.dat")
        const zfactor_filename = joinpath(res_directory, "zfactor_$(dim)d.dat")
    else
        const meff_filename = joinpath(res_directory, "meff_$(dim)d_spin$(spin).dat")
        const zfactor_filename = joinpath(res_directory, "zfactor_$(dim)d_spin$(spin).dat")
    end
    const sigma_z_filename = joinpath(data_directory, "data$(dim)d_Z_$(diagGenerate).jld2")
    const sigma_k_filename = joinpath(data_directory, "data$(dim)d_K_$(diagGenerate).jld2")
end
