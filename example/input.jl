dim = 2
rs = [1.0,]
# mass2 = [1.0, 2.0, 3.0]
mass2 = [0.05,]
Fs = [-0.0,]
beta = [25.0]
order = [3,]
neval = 1e6
# neval = 1e8
isDynamic = false
isFock = false

diagGenerate = :GV
# diagGenerate = :Parquet
isLayered2D = true

spin = 2
spinPolarPara = 2 / spin - 1 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]

if isLayered2D
    const parafilename = "para_wn_1minus0_layered2d.csv"
else
    const parafilename = "para_wn_1minus0.csv"
end