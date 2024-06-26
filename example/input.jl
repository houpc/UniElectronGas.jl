# Use finalized lambda scans to determine mass2 for maximum orders N = 4, 5, 6
include("lambda_scans.jl")
dim = 2      # dimension of the problem
rs = [2.0]
order = [4]  # maximum diagram order of the run
mass2 = rs_to_lambdas[dim][order[1]][rs[1]]

Fs = [-0.0]        # Fermi liquid parameter with zero angular momentum
beta = [50.0]      # inverse temperature beta = β*E_F 
neval = 2e8       # number of Monte Carlo samples
isDynamic = false  # whether to use effective field theory with dynamic screening or not 
isFock = false     # whether to use Fock renormalization or not

diagGenerate = :GV   # :GV or :Parquet, algorithm to generate diagrams
isLayered2D = true  # whether to use layered 2D system or not

spin = 2    # 2 for unpolarized, 1 for polarized
# spin = 1    # 2 for unpolarized, 1 for polarized
# spinPolarPara = 2 / spin - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
ispolarized = spin < 2

println("rs = $rs, mass2 = $mass2, order = $order, neval = $neval")

# Build file base names
basenames = [
    "para_wn_1minus0",
    "data$(dim)d_Z",
    "data$(dim)d_K",
    "data$(dim)d_dk",
    "meff_$(dim)d",
    "inverse_meff_$(dim)d",
    "zfactor_$(dim)d",
    "inverse_zfactor_$(dim)d",
    "chemical_potential_$(dim)d",
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
sigmadk_basename,
meff_basename,
inverse_meff_basename,
zfactor_basename,
inverse_zfactor_basename,
chemical_potential_basename,
dispersion_ratio_basename,
inverse_dispersion_ratio_basename = basenames

# Directory paths
data_directory = joinpath(@__DIR__, "sigma/data$(dim)d")
res_directory = joinpath(@__DIR__, "sigma")
para_directory = ""  # src directory

# File paths
const parafilename = joinpath(para_directory, para_basename * ".csv")
const sigma_z_filename = joinpath(data_directory, sigma_z_basename * ".jld2")
const sigma_k_filename = joinpath(data_directory, sigma_k_basename * ".jld2")
const sigmadk_filename = joinpath(data_directory, sigmadk_basename * ".jld2")
const meff_filename = joinpath(res_directory, meff_basename * ".dat")
const inverse_meff_filename = joinpath(res_directory, inverse_meff_basename * ".dat")
const zfactor_filename = joinpath(res_directory, zfactor_basename * ".dat")
const zinv_filename = joinpath(res_directory, inverse_zfactor_basename * ".dat")
const chemical_potential_filename = joinpath(res_directory, chemical_potential_basename * ".dat")
const dispersion_ratio_filename = joinpath(res_directory, dispersion_ratio_basename * ".dat")
const inverse_dispersion_ratio_filename = joinpath(res_directory, inverse_dispersion_ratio_basename * ".dat")
