import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import numpy as np


files = [
    "euler_a.csv", "euler_b.csv", "euler_c.csv",
    "e_cromer_a.csv", "e_cromer_b.csv", "e_cromer_c.csv",
    "midpoint_a.csv", "midpoint_b.csv", "midpoint_c.csv",
    "vscd_a.csv", "vscd_b.csv", "vscd_c.csv", "vscd_s.csv",
    "rk4_a.csv", "rk4_b.csv", "rk4_c.csv"
]

files_a = [
    "euler_a.csv",
    "e_cromer_a.csv",
    "midpoint_a.csv",
    "vscd_a.csv",
]

files_h = [
    "euler_h.csv", "e_cromer_h.csv",
    "midpoint_h.csv", "vscd_h.csv",
    "rk4_h.csv"
]

files_ac = [
    "ac_euler.csv", "ac_ecromer.csv",
    "ac_midpoint.csv", "ac_vscd.csv",
    "ac_rk4.csv"
]

files_ah = [
    "ah_euler.csv", "ah_ecromer.csv",
    "ah_midpoint.csv", "ah_vscd.csv",
    "ah_rk4.csv"
]

files_lle = [
    "lle_euler_a.csv", "lle_ecromer_a.csv",
    "lle_midpoint_a.csv", "lle_vscd_a.csv",
    "lle_rk4_a.csv"
]

def plot_bif_diagram(filename):
    data = pd.read_csv(filename)

    # убираем расширение
    name = os.path.splitext(filename)[0]
    # метод и параметр
    method, parameter = name.rsplit("_", 1)

    plt.figure(figsize=(12, 8))

    plt.scatter(
        data["parameter"],
        data[" peak"],
        s=1
    )

    plt.xlabel(f"Parameter {parameter}")
    plt.ylabel("Peak x")
    plt.title(f"Bifurcation diagram ({method}, parameter {parameter})")

    plt.show()


def plot_cluster_map(filename):
    data = pd.read_csv(filename)
    data.columns = data.columns.str.strip()

    x_name = data.columns[0]
    y_name = data.columns[1]
    z_name = data.columns[2]

    processed_data = []
    groups = data.groupby([x_name, y_name])

    for (x_val, y_val), group in groups:
        peaks = group[z_name].to_numpy().reshape(-1, 1)

        if len(peaks) == 0:
            count = 0
        else:
            db = DBSCAN(eps=0.01, min_samples=1).fit(peaks)
            count = len(set(db.labels_))

            if count > 15:
                count = 15

        processed_data.append({
            x_name: x_val,
            y_name: y_val,
            "period": count
        })

    result_df = pd.DataFrame(processed_data)

    plt.figure(figsize=(10, 8))

    scatter = plt.scatter(
        result_df[x_name],
        result_df[y_name],
        c=result_df["period"],
        cmap="turbo",
        s=7,
        marker='s'
    )

    cbar = plt.colorbar(scatter)
    cbar.set_label("Period (Number of unique peaks)")

    plt.xlabel(f"Parameter {x_name}")
    plt.ylabel(f"Parameter {y_name}")

    name = os.path.splitext(filename)[0]
    method = name.split("_", 1)[1] if "_" in name else name

    plt.title(f"2D Stability Map ({method})\nColor represents periodicity")

    plt.show()

def compare_bifurcations(file_method, file_reference):
    data_method = pd.read_csv(file_method)
    data_ref = pd.read_csv(file_reference)

    method = os.path.splitext(file_method)[0].split("_")[0]
    parameter = os.path.splitext(file_method)[0].split("_")[-1]

    plt.figure(figsize=(12, 8))

    plt.scatter(
        data_ref["parameter"],
        data_ref[" peak"],
        s=1,
        c="black",
        alpha=0.5,
        label="RK4 (reference)"
    )

    plt.scatter(
        data_method["parameter"],
        data_method[" peak"],
        s=1,
        c="red",
        alpha=0.5,
        label=method
    )

    plt.xlabel(f"Parameter {parameter}")
    plt.ylabel("Peak x")
    plt.title(f"{method} vs RK4")

    plt.legend()
    plt.show()


def plot_lle(filename):
    data = pd.read_csv(filename)

    plt.figure(figsize=(12, 6))

    plt.plot(
        data["parameter"],
        data["lambda"],
        linewidth=1.5
    )

    plt.axhline(0, color="red", linestyle="--", linewidth=1)

    plt.xlabel("Parameter")
    plt.ylabel("Largest Lyapunov Exponent")
    plt.title(f"Largest Lyapunov Exponent ({filename})")

    plt.grid(True)

    plt.show()

for file in files_lle:
    #plot_cluster_map(file)
    #compare_bifurcations(file, "rk4_a.csv")
    plot_lle(file)
