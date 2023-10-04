using ElectronLiquid
using ElectronLiquid.FeynmanDiagram

para = ParaMC(rs=1.0, beta=25, order=3, dim=3, isDynamic=false, isFock=false)
partition = UEG.partition(3)

diagram = Ver4.diagram(para, partition; channel=[PHr, PHEr, PPr], filter=[NoHartree, Proper])