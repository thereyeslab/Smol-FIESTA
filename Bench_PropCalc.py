#!/usr/bin/env python3
"""
compute_spot_metrics_csv.py

Under every "Analysis2" folder, find "*Cell_1_spots.csv",
compute in-Python the H→Q metrics (no Excel), and write out a new CSV
prefixed with "A" containing the original data plus the new columns.
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd

def compute_metrics(df):
    # G is column 7 (zero-based index 6)
    G = df.iloc[:, 6]

    # H = IF(G==0.001,1,0)
    df['H'] = (G == 0.001).astype(int)
    # I = IF(G==0.0011,1,0)
    df['I'] = (G == 0.0011).astype(int)

    # J = IF(AND(H==1, prev H==0),1,0)
    h = df['H'].to_numpy()
    prev_h = np.concatenate(([0], h[:-1]))
    df['J'] = ((h == 1) & (prev_h == 0)).astype(int)

    # K = IF(AND(I==1, prev I==0),1,0)
    i = df['I'].to_numpy()
    prev_i = np.concatenate(([0], i[:-1]))
    df['K'] = ((i == 1) & (prev_i == 0)).astype(int)

    # L = SUM(H)
    sum_h = int(df['H'].sum())
    df['L'] = sum_h

    # M = SUM(I)
    sum_i = int(df['I'].sum())
    df['M'] = sum_i

    # N = M/(M+L)
    df['N'] = sum_i / (sum_i + sum_h) if (sum_i + sum_h) else np.nan

    # O = SUM(J)
    sum_j = int(df['J'].sum())
    df['O'] = sum_j

    # P = SUM(K)
    sum_k = int(df['K'].sum())
    df['P'] = sum_k

    # Q = P/(P+O)
    df['Q'] = sum_k / (sum_k + sum_j) if (sum_k + sum_j) else np.nan

    return df

def process_csv(csv_path):
    print(f"⏳ Processing {csv_path}")
    df = pd.read_csv(csv_path)
    df = compute_metrics(df)

    # build output filename prefixed with "A"
    d, fname = os.path.split(csv_path)
    name, _ = os.path.splitext(fname)
    out_name = f"computed_{name}.csv"
    out_path = os.path.join(d, out_name)

    df.to_csv(out_path, index=False)
    print(f"✅ Saved: {out_path}")

def main(root_dir):
    # only look under Analysis2 folders now
    pattern = os.path.join(root_dir, '**', 'Analysis2', '*Cell_1_spots.csv')
    for csv_path in glob.iglob(pattern, recursive=True):
        process_csv(csv_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Compute H–Q spot metrics in Python and write new CSVs prefixed 'A', only under Analysis2 folders"
    )
    parser.add_argument('root_dir', help="Top‐level folder to search under")
    args = parser.parse_args()
    main(args.root_dir)
