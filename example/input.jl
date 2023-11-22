dim = 3 # dimension of the problem

### rs = 0.5 ###
# rs = [0.5]
# mass2 = [3.5]  # screening parameter

### rs = 1 ###
rs = [1.0]
mass2 = [1.75]  # screening parameter

### rs = 2 ###
# rs = [2.0]
# mass2 = [2.0]  # screening parameter

### rs = 3 ###
# rs = [3.0]
# mass2 = [1.5]  # screening parameter

### rs = 4 ###
# rs = [4.0]
# mass2 = [1.125]  # screening parameter

### rs = 5 ###
# rs = [5.0]
# mass2 = [1.125]  # screening parameter

### rs = 6 ###
# rs = [6.0]
# mass2 = [0.75]  # screening parameter

Fs = [-0.0]    # Fermi liquid parameter with zero angular momentum
beta = [40.0]   # inverse temperature beta = β*E_F 
order = [5]    # order of diagrams
neval = 1e11    # number of Monte Carlo samples
isDynamic = false # whether to use effective field theory with dynamic screening or not 
isFock = false # whether to use Fock renormalization or not

diagGenerate = :GV # :GV or :Parquet, algorithm to generate diagrams
isLayered2D = false # whether to use layered 2D system or not

spin = 2    # 2 for unpolarized, 1 for polarized
spinPolarPara = 2 / spin - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
ispolarized = spinPolarPara != 0.0

# Build file base names
basenames = [
    "para_wn_1minus0",
    "data$(dim)d_Z",
    "data$(dim)d_K",
    "meff_$(dim)d",
    "zfactor_$(dim)d",
    "inverse_zfactor_$(dim)d",
    "dispersion_ratio_$(dim)d",
    "inverse_dispersion_ratio_$(dim)d",
]
for i in eachindex(basenames)
    if spin != 2
        basenames[i] *= "_spin$(spin)"
    end
    if ispolarized
        basenames[i] *= "_polarized"
    end
    if isLayered2D
        @assert dim == 2 "Layered 2D mode is only available for dim = 2!"
        basenames[i] *= "_layered2d"
    end
end
para_basename,
sigma_z_basename,
sigma_k_basename,
meff_basename,
zfactor_basename,
inverse_zfactor_basename,
dispersion_ratio_basename,
inverse_dispersion_ratio_basename = basenames

# Directory paths
data_directory = joinpath(@__DIR__, "sigma/data$(dim)d")
res_directory = joinpath(@__DIR__, "sigma")
para_directory = ""  # src directory

# File paths
const parafilename = joinpath(para_directory, para_basename * ".csv")
const simga_z_filename = joinpath(data_directory, sigma_z_basename * ".jld2")
const sigma_k_filename = joinpath(data_directory, sigma_k_basename * ".jld2")
const meff_filename = joinpath(res_directory, meff_basename * ".dat")
const zfactor_filename = joinpath(res_directory, zfactor_basename * ".dat")
const zinv_filename = joinpath(res_directory, inverse_zfactor_basename * ".dat")
const dispersion_ratio_filename = joinpath(res_directory, dispersion_ratio_basename * ".dat")
const inverse_dispersion_ratio_filename = joinpath(res_directory, inverse_dispersion_ratio_basename * ".dat")
