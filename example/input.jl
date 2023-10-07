dim = 2. # dimension of the problem
rs = [2.0,]
# mass2 = [1.0, 2.0, 3.0]
mass2 = [0.1,]  # screening parameter
Fs = [-0.0,]    # Fermi liquid parameter with zero angular momentum
beta = [25.0]   # inverse temperature beta = β*E_F 
order = [3,]    # order of diagrams
neval = 1e6    # number of Monte Carlo samples
# neval = 1e8
isDynamic = false # whether to use effective field theory with dynamic screening or not 
isFock = false # whether to use Fock renormalization or not

diagGenerate = :GV # :GV or :Parquet, algorithm to generate diagrams
# diagGenerate = :Parquet
isLayered2D = true # whether to use layered 2D system or not

spin = 2    # 2 for unpolarized, 1 for polarized
spinPolarPara = 2 / spin - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]

if isLayered2D
    const parafilename = "para_wn_1minus0_layered2d.csv"
else
    const parafilename = "para_wn_1minus0.csv"
end