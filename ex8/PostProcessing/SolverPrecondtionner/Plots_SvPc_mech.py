#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

df = pd.read_csv("aggregation_SvPc_mech_1e5.csv")

df["preconditioner"] = df["preconditioner"].fillna("Aucun")

df_meca = df[(df["physics"] == "Mechanics") & (df["mean_s"] > 0)].copy()

if not df_meca.empty:
    df_meca = df_meca[df_meca["calls"] == 1].copy()

    mumps_data = df[(df["physics"] == "Mechanics") & (df["solver"] == "MUMPSSolver")]
    if not mumps_data.empty:
        mumps_ref = mumps_data["stat_mean"].iloc[0]
        if mumps_ref != 0:
            df_meca["rel_err"] = abs(df_meca["stat_mean"] - mumps_ref) / abs(mumps_ref)
            df_meca = df_meca[df_meca["rel_err"] <= 1e-8].copy()
        else:
            df_meca = df_meca[df_meca["stat_mean"] == 0].copy()

if not df_meca.empty:
    pivot_df = df_meca.pivot_table(index="solver", columns="preconditioner", values="mean_s")
    
    palette = ['#4C72B0', '#55A868', '#C44E52', '#8172B3', '#CCB974', '#64B5CD']
    color_dict = {}
    idx = 0
    for col in pivot_df.columns:
        if col == "Aucun":
            color_dict[col] = "#B3B3B3"
        else:
            color_dict[col] = palette[idx % len(palette)]
            idx += 1
            
    fig, ax = plt.subplots(figsize=(10, 7))
    
    bar_width = 0.15
    max_time = df_meca["mean_s"].max()
    ax.set_ylim(0, max_time * 1.35)
    
    for i, solver in enumerate(pivot_df.index):
        valid_data = pivot_df.loc[solver].dropna()
        n_bars = len(valid_data)
        
        if n_bars == 0:
            continue
            
        if n_bars == 1:
            offsets = [0]
        else:
            offsets = np.linspace(-bar_width * (n_bars - 1) / 2, bar_width * (n_bars - 1) / 2, n_bars)
            
        for offset, (precond, val) in zip(offsets, valid_data.items()):
            ax.bar(i + offset, val, width=bar_width, color=color_dict[precond], edgecolor="black", zorder=3)
            ax.text(i + offset, val + max_time * 0.02, f"{val:.1f}s", ha='center', va='bottom', rotation=90, fontsize=10)
            
    ax.set_title("Solving time : Mechanics (~1 million DOFs)", fontsize=16, weight="bold")
    ax.set_xlabel("Solver", fontsize=12, weight="bold")
    ax.set_ylabel("Time (s)", fontsize=12, weight="bold")
    
    ax.set_xticks(range(len(pivot_df.index)))
    ax.set_xticklabels(pivot_df.index, rotation=35, ha="right")
    ax.grid(axis="y", alpha=0.3, zorder=0)
    
    legend_elements = [Patch(facecolor=c, edgecolor='black', label=p) for p, c in color_dict.items()]
    ax.legend(handles=legend_elements, title="Preconditionner", fontsize=10, title_fontsize=11, loc="upper right")
    
    plt.tight_layout()
    plt.savefig("Profiling_SvPc_Mechanics_5e4.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    print("Graph has been generated : Profiling_SvPc_Mechanics_5e4.png")
else:
    print("No valid 'Mechanics' data was found in the CSV after filtering.")