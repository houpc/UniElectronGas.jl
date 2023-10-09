dim = 2 # dimension of the problem
rs = [1.0,]
mass2 = [1.0, 1.25, 1.5, 1.75, 2.0]  # screening parameter
Fs = [-0.0,]    # Fermi liquid parameter with zero angular momentum
beta = [40.0]   # inverse temperature beta = β*E_F 
order = [5,]    # order of diagrams
neval = 1e10    # number of Monte Carlo samples
isDynamic = false # whether to use effective field theory with dynamic screening or not 
isFock = false # whether to use Fock renormalization or not

diagGenerate = :GV # :GV or :Parquet, algorithm to generate diagrams
# diagGenerate = :Parquet
isLayered2D = false # whether to use layered 2D system or not

spin = 1    # 2 for unpolarized, 1 for polarized
spinPolarPara = 2 / spin - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
ispolarized = spinPolarPara != 0.0

data_directory = joinpath(@__DIR__, "sigma/data$(dim)d")

if isLayered2D
    const parafilename = "para_wn_1minus0_layered2d.csv"
    const simga_z_filename = joinpath(data_directory, "data$(dim)d_Z_$(diagGenerate)_layered.jld2")
    const sigma_k_filename = joinpath(data_directory, "data$(dim)d_K_$(diagGenerate)_layered.jld2")
    const zfactor_filename = spin == 2 ? "zfactor_layered2d.dat" : "zfactor_layered2d_spin$(spin).dat"
elseif ispolarized
    const parafilename = "para_wn_1minus0_GV_spin_polarized.csv"
    const sigma_z_filename = joinpath(data_directory, "data$(dim)d_Z_GV_spin_polarized.jld2")
    const sigma_k_filename = joinpath(data_directory, "data$(dim)d_K_GV_spin_polarized.jld2")
    const meff_filename = spin == 2 ? "meff_$(dim)d_GV_spin_polarized.dat" : "meff_$(dim)d_spin$(spin)_GV_spin_polarized.dat"
    const zfactor_filename = spin == 2 ? "zfactor_$(dim)d_GV_spin_polarized.dat" : "zfactor_$(dim)d_spin$(spin)_GV_spin_polarized.dat"
else
    const parafilename = "para_wn_1minus0.csv"
    const sigma_z_filename = joinpath(data_directory, "data$(dim)d_Z_$(diagGenerate).jld2")
    const sigma_k_filename = joinpath(data_directory, "data$(dim)d_K_$(diagGenerate).jld2")
    const meff_filename = spin == 2 ? "meff_$(dim)d.dat" : "meff_$(dim)d_spin$(spin).dat"
    const zfactor_filename = spin == 2 ? "zfactor_$(dim)d.dat" : "zfactor_$(dim)d_spin$(spin).dat"
end
