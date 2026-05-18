#!/usr/bin/env python3
import os
import sys
import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def get_particle_label(tree_name):
    particle_labels = {
        "X3872": "X(3872)",
        "Psi2S": "psi(2S)",
        "ntKp": "B+",
        "ntKstar": "B0",
        "ntphi": "Bs",
    }
    return particle_labels.get(tree_name, tree_name)
    

def resolve_paths(tree_name, system_name):
    base_dir = "/eos/user/h/hmarques/Analysis_CODES"
    x3872_sample_dir = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872"

    if tree_name in ("X3872", "Psi2S"):
        sample_dirs = {
            "ppRef": ("ppRef24", "ppRef"),
            "ppRef24": ("ppRef24", "ppRef"),
            "PbPb23": ("PbPb23", "PbPb23"),
            "PbPb24": ("PbPb24", "PbPb24"),
        }
        if system_name not in sample_dirs:
            raise ValueError(f"Unknown X3872/Psi2S system '{system_name}'")

        sample_subdir, sample_tag = sample_dirs[system_name]
        data_tree = "ntmix"
        mc_tree = "ntmix_PSI2S" if tree_name == "Psi2S" else "ntmix_X3872"
        path_to_data = f"{x3872_sample_dir}/{sample_subdir}/flat_ntmix_{sample_tag}_DATA.root"
        path_to_mc = (
            f"{x3872_sample_dir}/{sample_subdir}/flat_ntmix_{sample_tag}_MC_PSI2S.root"
            if tree_name == "Psi2S"
            else f"{x3872_sample_dir}/{sample_subdir}/flat_ntmix_{sample_tag}_MC_X3872.root"
        )
    else:
        data_tree = tree_name
        mc_tree = tree_name
        path_to_data = f"{base_dir}/flatER/Bmeson/flat_{tree_name}_{system_name}_DATA.root"
        path_to_mc = f"{base_dir}/flatER/Bmeson/flat_{tree_name}_{system_name}_MC.root"

    return path_to_data, path_to_mc, data_tree, mc_tree


def plot_correlation_dataMC(tree="X3872", system_name="ppRef"):
    ROOT.gROOT.SetBatch(True)

    variables = [
        "Btrk1Pt",
        "Btrk1dR",
        "Btrk2dR",
        "BtrkPtimb",
        "BtktkvProb",
        "Bchi2Prob",
        "Bmass",
        "BQvalue",
        "Bpt",
    ]

    path_to_data, path_to_mc, data_tree_name, mc_tree_name = resolve_paths(tree, system_name)

    selection = "1"
    sideband = "1"
    if data_tree_name == "ntmix":
        sideband = "(Bmass > 3.95 || (Bmass > 3.75 && Bmass < 3.8))"
    else:
        sideband = "(Bmass > 5.5)"

    data_cut = f"({selection}) && ({sideband})"
    data_file = ROOT.TFile.Open(path_to_data)
    mc_file = ROOT.TFile.Open(path_to_mc)
    data_tree = data_file.Get(data_tree_name)
    mc_tree = mc_file.Get(mc_tree_name)

    common_vars = []
    for var in variables:
        has_data = data_tree.GetListOfBranches().FindObject(var)
        has_mc = mc_tree.GetListOfBranches().FindObject(var)
        if has_data and has_mc:
            common_vars.append(var)

    print("CASE:", tree)
    print("systemNAME:", system_name)
    print("Selection:", selection)
    print("Sideband:", sideband)

    rdf_data = ROOT.RDataFrame(data_tree_name, path_to_data).Filter(data_cut)
    arr_data = rdf_data.AsNumpy(common_vars)
    n_data = len(arr_data[common_vars[0]])
    mat_data = np.column_stack([np.asarray(arr_data[var], dtype=float) for var in common_vars])
    corr_data = np.corrcoef(mat_data, rowvar=False)
    if data_tree_name == "ntmix":
        out_dir = "./ntmix"
        data_output_tag = "ntmix"
    else:
        out_dir = f"./{tree}"
        data_output_tag = tree
    os.makedirs(out_dir, exist_ok=True)

    draw_matrix(
        corr_data,
        common_vars,
        f"data sideband {system_name}",
        f"{out_dir}/data_sideband_corr_{system_name}_{data_output_tag}.pdf",
    )

    rdf_mc = ROOT.RDataFrame(mc_tree_name, path_to_mc).Filter(selection)
    arr_mc = rdf_mc.AsNumpy(common_vars)
    n_mc = len(arr_mc[common_vars[0]])
    print("n_data:", n_data)
    print("n_mc:", n_mc)
    mat_mc = np.column_stack([np.asarray(arr_mc[var], dtype=float) for var in common_vars])
    corr_mc = np.corrcoef(mat_mc, rowvar=False)
    draw_matrix(
        corr_mc,
        common_vars,
        f"{get_particle_label(tree)} signal MC {system_name}",
        f"{out_dir}/mc_corr_{system_name}_{tree}.pdf",
    )


def draw_matrix(matrix, labels, title, out_name, vmax=1.0):
    fig_size = max(8, int(0.7 * len(labels)))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    ax.set_title(title, fontsize=20, fontweight="bold")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=8)

    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_name)
    plt.close(fig)
    print(f"Saved correlation matrix to {out_name}")


if __name__ == "__main__":
    tree = "X3872"
    system = "ppRef"
    if len(sys.argv) > 1:
        tree = sys.argv[1]
    if len(sys.argv) > 2:
        system = sys.argv[2]
    plot_correlation_dataMC(tree, system)


## python3 correlation_dataMC.py X3872 ppRef
## python3 correlation_dataMC.py Psi2S ppRef
## python3 correlation_dataMC.py ntphi ppRef
