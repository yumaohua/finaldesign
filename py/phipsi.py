#!/usr/bin/env python3
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_delta(csv_path, out_dir="pic"):
    # 读入 CSV
    df = pd.read_csv(csv_path)
    # 生成画布
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(df['idx'], df['delta_phi_re'], label='Δφ real')
    ax.plot(df['idx'], df['delta_phi_im'], label='Δφ imag', linestyle='--')
    ax.plot(df['idx'], df['delta_psi_re'], label='Δψ real')
    ax.plot(df['idx'], df['delta_psi_im'], label='Δψ imag', linestyle='--')
    ax.set_xlabel('Element index $l$')
    ax.set_ylabel('Δ value')
    ax.set_title(os.path.basename(csv_path))
    ax.legend(loc='best')
    ax.grid(True)
    fig.tight_layout()

    # 确保输出目录存在
    os.makedirs(out_dir, exist_ok=True)
    # 保存到 pic/ 目录，文件名与 CSV 同名但后缀 .png
    png_name = os.path.splitext(os.path.basename(csv_path))[0] + '.png'
    out_path = os.path.join(out_dir, png_name)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")

def main():
    # 从 data/ 目录读取 delta_*.csv
    csv_files = sorted(glob.glob('./data/delta_*.csv'))
    if not csv_files:
        print("No CSV files found in data/")
        return
    for csv in csv_files:
        plot_delta(csv, out_dir="./pic")

if __name__ == '__main__':
    main()
