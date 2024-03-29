# dim = 2 # dimension of the problem
dim = 3 # dimension of the problem
# rs = [0.5]
# rs = [5.0]
rs = [1.0]
# mass2 = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,]  # screening parameter
# mass2 = [1.125]  # screening parameter
mass2 = [1.0,]  # screening parameter
Fs = [-0.0,]    # Fermi liquid parameter with zero angular momentum
# beta = [10.0, 20.0, 40.0, 80.0, 160.0]   # inverse temperature beta = β*E_F 
# beta = [5.0, 10.0, 20.0, 40.0, 80.0, 160.0]   # inverse temperature beta = β*E_F 
# beta = [80.0,]   # inverse temperature beta = β*E_F 
beta = [25.0,]   # inverse temperature beta = β*E_F 
order = [4,]    # order of diagrams
# order = [5,]    # order of diagrams
# neval = 5e7    # number of Monte Carlo samples
# order = [6,]    # order of diagrams
neval = 1e7    # number of Monte Carlo samples
isDynamic = false # whether to use effective field theory with dynamic screening or not 
isFock = false # whether to use Fock renormalization or not

# diagGenerate = :GV # :GV or :Parquet, algorithm to generate diagrams
diagGenerate = :Parquet
# isLayered2D = true # whether to use layered 2D system or not
isLayered2D = false # whether to use layered 2D system or not

spin = 2    # 2 for unpolarized, 1 for polarized
# spin = 1    # 2 for unpolarized, 1 for polarized
spinPolarPara = 2 / spin - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]

if isLayered2D
    const parafilename = "para_wn_1minus0_layered2d.csv"
else
    # const parafilename = "para_wn_1minus0.csv"
    const parafilename = "para_wn_1minus0_3d.csv"
end
