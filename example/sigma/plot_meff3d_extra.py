#!/usr/bin/python3
import scienceplots
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib as mat
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

cdict = {
    "blue": "#0077BB",
    "cyan": "#33BBEE",
    "teal": "#009988",
    "orange": "#EE7733",
    "red": "#CC3311",
    "magenta": "#EE3377",
    "grey": "#BBBBBB",
}
# plt.switch_backend("TkAgg")
plt.style.use(["science", "std-colors"])
# mat.rcParams.update({'font.size': 28})
# mat.rcParams.update({"font.size": 16})
mat.rcParams["font.size"] = 16
mat.rcParams["mathtext.fontset"] = "cm"
mat.rcParams["font.family"] = "Times New Roman"
# size = 36
# sizein = 12

points = ['o', '^', 's', 'v', 'p', '<', 'h', '>']
# Colormap = ["blue", "red", "orange", "violet", "green", "black", "darkblue"]
# Color_beta = {0.125: "red", 0.25: "orange", 0.5: "violet", 1.0: "green", 2.0: "blue", 2.5:"darkblue", 16.0: "black"}
# Colormap = [
#     cdict["blue"],
#     cdict["red"],
#     cdict["orange"],
#     cdict["magenta"],
#     cdict["cyan"],
#     cdict["teal"],
#     "grey",
#     "black",
# ]
Colormap = [
    cdict["blue"],
    cdict["red"],
    cdict["orange"],
    cdict["magenta"],
    "grey",
    cdict["teal"],
    cdict["cyan"],
    "black",
]
# pt_beta = {0.125: 'o', 0.25: '^', 0.5: 'h', 1.0: 'v', 2.0: 'd', 2.5:'p', 16.0: 's'}
pts = ['s', '^', 'v', 'p', 's', '<', 'h', 'o']
# fig, axes = plt.subplots(1,2, sharex='col')#, sharey=True, gridspec_kw={'wspace': 0})
# fig, ax1 = plt.subplots(figsize=(6, 4.6))
# fig, ax1 = plt.subplots(figsize=(6, 4))
fig, ax1 = plt.subplots(figsize=(5, 5))
# fig.set_size_inches(13,10)

# m_VDMC = [1.0, 0.95861, 0.9511, 0.9516, 0.9601, 0.968, 0.9782]
# m_VDMC_err = [0, 0.00089, 0.0021, 0.0031, 0.0022, 0.0021, 0.0025]
# rs_VDMC = [0, 0.5, 1, 2, 3, 4, 5]
m_VDMC = [1.0, 0.95893, 0.9514, 0.9516, 0.9597, 0.9692, 0.9771, 0.9842]
m_VDMC_err = [0, 0.00067, 0.0016, 0.0018, 0.0016, 0.0026, 0.0028, 0.0029]
rs_VDMC = [0, 0.5, 1, 2, 3, 4, 5, 6]

rs_VMC = [0, 1, 2, 4, 5, 10]
# m_BFVMC = [[1.00, 0.01], [0.98, 0.01], [1.00, 0.02], [1.09,0.03],[1.28, 0.03]]
m_BFVMC = [1.0, 1.00, 0.98, 1.00, 1.09, 1.28]
m_BFVMC_err = [0, 0.01, 0.01, 0.02, 0.03, 0.03]
# m_SJVMC = [[0.96,0.01], [0.94, 0.02], [0.94, 0.02], [1.02, 0.02], [1.13, 0.03]]
m_SJVMC = [1.0, 0.96, 0.94, 0.94, 1.02, 1.13]
m_SJVMC_err = [0, 0.01, 0.02, 0.02, 0.02, 0.03]

rs_DMC = [0, 1, 2, 3, 4, 5]
m_DMC = [1.0, 0.918, 0.879, 0.856, 0.842, 0.791]
m_DMC_err = [0, 0.006, 0.014, 0.014, 0.017, 0.01]

f_sigma_G0W0 = np.load(
    "results/finalized_meff_results/3d/rpa/meff_3d_sigma_G0W0.npz")
rs_G0W0 = f_sigma_G0W0["rslist"]
m_G0W0 = f_sigma_G0W0["mefflist"]

f_sigma_G0Wp = np.load(
    "results/finalized_meff_results/3d/simion_giuliani/meff_3d_sigma_G0Wp.npz")
# f_sigma_G0Wp = np.load(
#     "results/finalized_meff_results/3d/const/meff_3d_sigma_G0Wp.npz")
rs_G0Wp = f_sigma_G0Wp["rslist"]
m_G0Wp = f_sigma_G0Wp["mefflist"]
# TODO: Fix these outliers
mask = [2, 3]
rs_G0Wp = np.delete(rs_G0Wp, mask)
m_G0Wp = np.delete(m_G0Wp, mask)

f_sigma_G0Wpm = np.load(
    "results/finalized_meff_results/3d/simion_giuliani/meff_3d_sigma_G0Wpm.npz")
# f_sigma_G0Wpm = np.load(
#     "results/finalized_meff_results/3d/const/meff_3d_sigma_G0Wpm.npz")
rs_G0Wpm = f_sigma_G0Wpm["rslist"]
m_G0Wpm = f_sigma_G0Wpm["mefflist"]

f_tree_level_G0W0 = np.load(
    "results/finalized_meff_results/3d/rpa/meff_3d_tree_level_G0W0.npz")
rs_tree_level_G0W0 = f_tree_level_G0W0["rslist"]
m_tree_level_G0W0 = f_tree_level_G0W0["mefflist"]

f_tree_level_G0Wp = np.load(
    "results/finalized_meff_results/3d/const/meff_3d_tree_level_G0Wp.npz")
rs_tree_level_G0Wp = f_tree_level_G0Wp["rslist"]
m_tree_level_G0Wp = f_tree_level_G0Wp["mefflist"]

rs_QSGW = [2, 3, 4, 5]
m_QSGW = [1.0072289156626506, 1.0188253012048194, 1.0335843373493976, 1.0562500000000001]

# rs_RPA = [0, 1, 2, 3, 4, 5, 6]
# # m_RPA = [[0.97, 0], [0.992,0], [1.016,0] [1.039,0], [1.059, 0]]
# # m_Gp = [[0.952, 0], [0.951, 0], [0.956, 0], [0.962, 0], [0.968, 0]]
# # m_Gpm = [[0.957, 0], [0.966, 0], [0.983, 0], [1.005, 0], [1.028,0]]
# m_G0W0 = [1.0, 0.970, 0.992, 1.016, 1.039, 1.059, 1.078]
# m_Gp = [1.0, 0.952, 0.951, 0.956, 0.962, 0.968, 0.973]
# m_Gpm = [1.0, 0.957, 0.966, 0.983, 1.005, 1.028, 1.055]


def errorbar_mvsrs(rs, meff_data, merr, idx, label, ax=ax1, zorder=None, capsize=4):
    if zorder is not None:
        handle = ax.errorbar(rs, meff_data, merr, fmt=pts[idx], capthick=1, capsize=capsize,
                             ms=5, color=Colormap[idx], label=label, zorder=zorder)
    else:
        handle = ax.errorbar(rs, meff_data, merr, fmt=pts[idx], capthick=1, capsize=capsize,
                             ms=5, color=Colormap[idx], label=label)
    return handle


def plot_mvsrs(rs, meff_data, idx, label, ls='-', ax=ax1, rs_HDL=None, meff_HDL=None, zorder=None):
    # Add data in the high-density limit to the fit, if provided
    if (rs_HDL is not None) and (meff_HDL is not None):
        rslist = pd.unique(np.concatenate([rs, rs_HDL]))
        mdatalist = pd.unique(np.concatenate([meff_data, meff_HDL]))
        P = np.argsort(rslist)
        rs = rslist[P]
        meff_data = mdatalist[P]
    print(rs)
    print(meff_data)

    # mfitfunc = interpolate.PchipInterpolator(rs, meff_data)
    mfitfunc = interpolate.Akima1DInterpolator(rs, meff_data)
    # xgrid = np.arange(0, 6.2, 0.02)
    xgrid = np.arange(0, 6.2, 0.02)
    # ax.plot(rs, meff_data, 'o', ms=10, color=Colormap[idx])
    # ax.plot(xgrid, mfitfunc(xgrid), ls=ls, lw=2, color=Colormap[idx], label=label)
    if zorder is not None:
        handle, = ax.plot(xgrid, mfitfunc(xgrid), ls=ls,
                          color=Colormap[idx], label=label, zorder=zorder)
    else:
        handle, = ax.plot(xgrid, mfitfunc(xgrid), ls=ls,
                          color=Colormap[idx], label=label)

    yfit = np.ma.masked_invalid(mfitfunc(xgrid))
    # print(yfit)
    print("Turning point: rs = ", xgrid[np.argmin(yfit)])
    print("Effective mass ratio at turning point: ", np.min(yfit))
    return handle

# Taylor series for m* / m in the high-density limit to leading order in rs in 3D
# (c.f. Giuliani & Vignale, Quantum Theory of the Electron Liquid, 2008, p. 500)
def high_density_limit(rs):
    alpha = (4.0 / (9.0 * np.pi))**(1.0 / 3.0)
    return 1 + alpha * rs * np.log(rs) / (2.0 * np.pi)


# High-density limit
rs_HDL_plot = np.linspace(1e-5, 0.35, num=101)
meff_HDL_plot = np.array([high_density_limit(rs) for rs in rs_HDL_plot])

# Use exact expression in the high-density limit for RPA(+FL) and RPT fits
cutoff_HDL = 0.1
rs_HDL = rs_HDL_plot[rs_HDL_plot <= cutoff_HDL]
meff_HDL = meff_HDL_plot[rs_HDL_plot <= cutoff_HDL]
print(rs_HDL)
print(meff_HDL)

# High-density limit
# plot_mvsrs(rs_HDL_plot, meff_HDL_plot, 7, "HDL", "--", zorder=1000)

# For paper
# reflabels = [" [12]", " [13]"]

# For talk
reflabels = ["$^*$", "$^\dagger$", "$^\ddagger$"]

# Digitized data from Kutepov 3DUEG LQSGW paper
# handle1 = ax1.scatter(rs_QSGW, m_QSGW, 20, color=Colormap[0], label=rf"QSGW{reflabels[0]}", marker='s', zorder=30)

# m_QMC = [m_DMC, m_SJVMC]
# m_QMC_errs = [m_DMC_err, m_SJVMC_err]
# labels = ["DMC (2021)", "VMC (2023)"]
# idx = 0
# for (mdat, merr, label) in zip(m_QMC, m_QMC_errs, labels):
# errorbar_mvsrs(rs_VMC, mdat, merr, idx, label)
# idx = idx + 1
# handle1 = errorbar_mvsrs(rs_DMC, m_DMC, m_DMC_err, 0,
#                          rf"DMC{reflabels[0]}", zorder=10)
handle2 = errorbar_mvsrs(rs_VMC, m_SJVMC, m_SJVMC_err,
                         1, rf"VMC{reflabels[1]}", zorder=20)
#  1, rf"VMC{reflabels[1]}", zorder=20)

idx = 2
l1_handles = []
# rsdats = [rs_G0W0, rs_tree_level_G0W0, rs_tree_level_G0Wp]
# mdats = [m_G0W0, m_tree_level_G0W0, m_tree_level_G0Wp]
# labels = [r"RPA", r"Tree-level $G_0 W_0$", r"Tree-level $G_0 W_{+}$"]
rsdats = [rs_tree_level_G0W0, rs_tree_level_G0Wp]
mdats = [m_tree_level_G0W0, m_tree_level_G0Wp]
labels = [r"Tree-level $G_0 W_0$", r"Tree-level $G_0 W_{+}$"]
for (rsdat, mdat, label) in zip(rsdats, mdats, labels):
    print("\nPlotting ", label)
    handle = plot_mvsrs(rsdat, mdat, idx, label, '--',
                        rs_HDL=rs_HDL, meff_HDL=meff_HDL)
    l1_handles.append(handle)
    idx = idx + 1

# l1_handles = []
# m_RPA = [m_G0W0]
# rs_RPA = [rs_G0Wp, rs_G0Wpm]
# m_RPA = [m_G0Wp, m_G0Wpm]
rs_RPA = [rs_G0W0, rs_G0Wp, rs_G0Wpm]
m_RPA = [m_G0W0, m_G0Wp, m_G0Wpm]
# merr = np.zeros(len(rs_RPA))
# labels = [r"RPA", rf"$G_0 W_+$", rf"$G_0 W_\pm$"]
labels = [r"RPA", rf"$G_0 W_+${reflabels[0]}", rf"$G_0 W_\pm${reflabels[0]}"]
# labels = [r"RPA"]
# labels = [r"$G_0 W_0$"]
# labels = [r"$G_0 W_0$", r"$G_+$ [15]", r"$G_+$\,\&\,$G_-$ [15]"]
for (rsdat, mdat, label) in zip(rs_RPA, m_RPA, labels):
    print("\nPlotting ", label)
    handle = plot_mvsrs(rsdat, mdat, idx, label, '--',
                        rs_HDL=rs_HDL, meff_HDL=meff_HDL)
    l1_handles.append(handle)
    idx = idx + 1

# Add points at rs = 0
# rs_HDL = np.insert(rs_HDL, 0, 0.0)
# meff_HDL = np.insert(meff_HDL, 0, 1.0)

handle3 = errorbar_mvsrs(rs_VDMC, m_VDMC, m_VDMC_err,
                         7, "This work", zorder=30)
# rs_VDMC.append(5.2)
# m_VDMC.append(0.98)
print("\nPlotting our data")
plot_mvsrs(rs_VDMC, m_VDMC, 7, "", '-', rs_HDL=rs_HDL, meff_HDL=meff_HDL)

ax1.set_xlabel(r"$r_s$")
ax1.set_ylabel(r"$m^*/m$")
ax1.set_xlim(0, 6.2)
ax1.set_ylim(0.875, 1.105)
# ax1.set_ylim(0.83, 1.14)
# ax1.set_ylim(0.765, 1.135)
# ax1.set_ylim(0.765, 1.145)
ax1.annotate(r"3D", xy=(0.875, 0.9), xycoords="axes fraction")

# Assemble legends
l2_handles = [handle2, handle3]
# l2_handles = [handle1, handle2, handle3]
top_legend = plt.legend(handles=l1_handles, loc="upper left", fontsize=14)
bottom_legend = plt.legend(handles=l2_handles, loc="lower left", fontsize=14)
ax1.add_artist(top_legend)
ax1.add_artist(bottom_legend)
ax1.set_xticks([0, 1, 2, 3, 4, 5, 6])
# ax1.set_yticks([0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1])
# ax1.set_yticks([0.95, 1.0, 1.05])

# Save the figure
# plt.savefig("figures/meff_3DUEG.pdf")

# For talk
# plt.savefig("figures/meff_3DUEG_talk.pdf")
plt.savefig("figures/meff_3DUEG_extra.pdf")
# plt.savefig("figures/meff_3DUEG_extra_const.pdf")
# plt.savefig("figures/meff_3DUEG_extra_simion_giuliani.pdf")
