#!/usr/bin/env python
"""
Runnable Script if run as __main__
Perform exponential mixture model fitting (with photobleaching correction)
to the “rebind‐strict‐boundtime.csv” file produced by your analysis pipeline.
The input file is constructed as:
    {csv_path}/{output_folder_name}/rebind‑strict‑boundtime.csv
The measurement column is assumed to be "Bound Time" (in seconds after scaling by frameRate).
An interval of 10 ms is assumed (so frameRate = 100) and the photobleaching rate is fixed.
The script fits both a single‐ and a double‐exponential model and saves the results
and corresponding plots to a PDF file in the output folder.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import minimize
from scipy.stats import expon
from sklearn.cluster import KMeans
import tomllib
import argparse

# Set number of threads to 1 for reproducibility
os.environ["OMP_NUM_THREADS"] = "1"


# ---------------------- Negative Log-Likelihood Functions ----------------------

def neg_log_likelihood_single(params, data):
    lambda1 = params[0]
    likelihoods = expon.pdf(data, scale=1 / lambda1)
    neg_log_lik = -np.sum(np.log(likelihoods))
    return neg_log_lik


def neg_log_likelihood_with_bleaching(params, data, lambda_bleach):
    lambda1, lambda2, p1 = params
    p2 = 1 - p1
    # In this model the photobleaching rate is assumed fixed externally.
    likelihoods = (p1 * expon.pdf(data, scale=1 / lambda1) +
                   p2 * expon.pdf(data, scale=1 / lambda2))
    neg_log_lik = -np.sum(np.log(likelihoods))
    return neg_log_lik


# ---------------------- Initialization ----------------------

def get_initial_parameters(data):
    data_reshaped = data.reshape(-1, 1)
    kmeans = KMeans(n_clusters=2, n_init=10, random_state=0).fit(data_reshaped)
    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_.flatten()
    lambda1_initial = 1 / cluster_centers[0]
    lambda2_initial = 1 / cluster_centers[1]
    p1_initial = np.sum(labels == 0) / len(data)
    return lambda1_initial, lambda2_initial, p1_initial


# ---------------------- Model Fitting Functions ----------------------

def run_single_exponential_model(data):
    initial_params = [1 / np.mean(data)]
    result = minimize(neg_log_likelihood_single, initial_params, args=(data,), bounds=[(1e-6, None)])
    lambda1 = result.x[0]
    neg_log_lik = neg_log_likelihood_single([lambda1], data)
    n = len(data)
    bic = 2 * neg_log_lik + np.log(n) * 1  # 1 parameter
    mean_time_single = 1 / lambda1
    return lambda1, mean_time_single, bic


def run_mixture_model(data_path, lambda_bleach, measurement, frameRate, max_iterations=100, tolerance=0.05):
    # Read the CSV file from the analysis pipeline
    data = pd.read_csv(data_path)[measurement].values / frameRate
    # (Here the division by frameRate converts the raw "Bound Time" into seconds)
    for i in range(max_iterations):
        lambda1_initial, lambda2_initial, p1_initial = get_initial_parameters(data)
        initial_params = [lambda1_initial, lambda2_initial, p1_initial]
        result = minimize(neg_log_likelihood_with_bleaching, initial_params,
                          args=(data, lambda_bleach),
                          bounds=[(1e-6, None), (1e-6, None), (0, 1)])
        lambda1, lambda2, p1 = result.x
        p2 = 1 - p1
        # Ensure lambda1 corresponds to the slower decay (i.e. larger mean)
        if lambda1 < lambda2:
            lambda1, lambda2 = lambda2, lambda1
            p1, p2 = p2, p1
        # Compute corrected mean times (excluding photobleaching effect)
        try:
            mean_time1_corrected = (1 / lambda1) * (1 / lambda_bleach) / ((1 / lambda_bleach) - (1 / lambda1))
            mean_time2_corrected = (1 / lambda2) * (1 / lambda_bleach) / ((1 / lambda_bleach) - (1 / lambda2))
        except ZeroDivisionError:
            mean_time1_corrected, mean_time2_corrected = np.inf, np.inf
        # Non-corrected mean times (including photobleaching effect)
        mean_time1_non_corrected = 1 / lambda1
        mean_time2_non_corrected = 1 / lambda2
        # (We iterate max_iterations times for now; you could add a convergence check if desired.)
    neg_log_lik = neg_log_likelihood_with_bleaching([lambda1, lambda2, p1], data, lambda_bleach)
    n = len(data)
    bic = 2 * neg_log_lik + np.log(n) * 3  # 3 parameters
    # Calculate weight fractions based on non-corrected means
    weight_fraction1 = (mean_time1_non_corrected * p1) / (mean_time1_non_corrected * p1 + mean_time2_non_corrected * p2)
    weight_fraction2 = (mean_time2_non_corrected * p2) / (mean_time1_non_corrected * p1 + mean_time2_non_corrected * p2)
    results = {
        "Estimated Mean Time 1 (corrected)": mean_time1_corrected,
        "Estimated Decay Rate 1 (corrected)": 1 / mean_time1_corrected if mean_time1_corrected != 0 else np.inf,
        "Estimated Mean Time 2 (corrected)": mean_time2_corrected,
        "Estimated Decay Rate 2 (corrected)": 1 / mean_time2_corrected if mean_time2_corrected != 0 else np.inf,
        "Estimated Mean Time 1 (non-corrected)": mean_time1_non_corrected,
        "Estimated Decay Rate 1 (non-corrected)": lambda1,
        "Estimated Mean Time 2 (non-corrected)": mean_time2_non_corrected,
        "Estimated Decay Rate 2 (non-corrected)": lambda2,
        "Fraction of Exponential Component 1": p1,
        "Fraction of Exponential Component 2": p2,
        "Weight Fraction 1": weight_fraction1,
        "Weight Fraction 2": weight_fraction2,
        "BIC Double Exponential": bic
    }
    return results, data, lambda1, lambda2, p1, p2


# ---------------------- Plotting Functions ----------------------

def plot_results(data, lambda1, lambda2, p1, p2, lambda1_single=None, mean_time_single=None, measurement="Time"):
    bins = np.linspace(0, data.max(), 50)
    hist, bins = np.histogram(data, bins=bins, density=True)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, density=True, alpha=0.6, color='grey', label='Data', edgecolor='black')
    if lambda1_single is not None and mean_time_single is not None:
        pdf_single = expon.pdf(bin_centers, scale=1 / lambda1_single)
        plt.plot(bin_centers, pdf_single, 'g--', linewidth=2,
                 label=f'Single Exp Fit (mean={mean_time_single:.3f}, λ={lambda1_single:.3f})')
    else:
        pdf1 = p1 * expon.pdf(bin_centers, scale=1 / lambda1)
        pdf2 = p2 * expon.pdf(bin_centers, scale=1 / lambda2)
        plt.fill_between(bin_centers, pdf2, pdf1 + pdf2, color='red', alpha=0.5,
                         label=f'Exp 1 (mean={1 / lambda1:.3f}, λ={lambda1:.3f})')
        plt.fill_between(bin_centers, 0, pdf2, color='blue', alpha=0.5,
                         label=f'Exp 2 (mean={1 / lambda2:.3f}, λ={lambda2:.3f})')
        plt.plot(bin_centers, pdf1 + pdf2, 'k-', linewidth=2, label='Total Double Exp Fit')
    plt.xlabel(measurement + " (seconds)")
    plt.ylabel("Density")
    plt.legend()
    plt.title("Histogram with Exponential Fits")
    plt.tight_layout()


def plot_results_with_trendlines(data, lambda1, lambda2, p1, p2, lambda1_single=None, mean_time_single=None,
                                 measurement="Time"):
    bins = np.linspace(0, data.max(), 50)
    hist, bins = np.histogram(data, bins=bins, density=True)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=bins, density=True, alpha=0.6, color='grey', label='Data', edgecolor='black')
    pdf1 = p1 * expon.pdf(bin_centers, scale=1 / lambda1)
    pdf2 = p2 * expon.pdf(bin_centers, scale=1 / lambda2)
    plt.plot(bin_centers, pdf1, 'r-', linewidth=2,
             label=f'Exp 1 Trendline (mean={1 / lambda1:.3f}, λ={lambda1:.3f})')
    plt.plot(bin_centers, pdf2, 'b-', linewidth=2,
             label=f'Exp 2 Trendline (mean={1 / lambda2:.3f}, λ={lambda2:.3f})')
    plt.xlabel(measurement + " (seconds)")
    plt.ylabel("Density")
    plt.legend()
    plt.title("Histogram with Double Exponential Trendlines")
    plt.tight_layout()


# ---------------------- Main ----------------------

def main(config_path: str = None):
    # Read configuration from TOML
    if not config_path:
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, "script-config.toml")
    with open(config_path, "rb") as config_file:
        configs = tomllib.load(config_file)

    # Get paths from config
    csv_path = configs["path"]["csv_path"]
    output_folder = configs["path"]["output_folder_name"]
    # Construct the data path for the rebind strict bound time file
    data_path = os.path.join(csv_path, output_folder, "rebind-strict-boundtime.csv")

    # Define fixed parameters (no new configs)
    measurement = "Bound Time"  # Column name in CSV
    interval_ms = 10  # in milliseconds
    frameRate = 1000 / interval_ms  # Conversion factor (e.g., 100 frames per second)
    lambda_bleach = 1 / 10000000  # Fixed photobleaching rate

    # Run the mixture model
    results, data, lambda1, lambda2, p1, p2 = run_mixture_model(data_path, lambda_bleach, measurement, frameRate)

    # Run the single exponential model
    lambda1_single, mean_time_single, bic_single = run_single_exponential_model(data)
    results["BIC Single Exponential"] = bic_single
    results["Estimated Mean Time (Single Exponential)"] = mean_time_single
    results["Estimated Decay Rate (Single Exponential)"] = lambda1_single

    # Calculate weight fraction for single exponential (trivial in this case)
    weight_fraction_single = mean_time_single / (mean_time_single * len(data))
    results["Weight Fraction (Single Exponential)"] = weight_fraction_single

    # Save results and plots to a single PDF in the same output folder
    pdf_path = os.path.join(csv_path, output_folder, "Exponential_fit_diffusion_B.pdf")
    with PdfPages(pdf_path) as pdf:
        # Save text results as a figure
        plt.figure(figsize=(8, 4))
        textstr = "\n".join([f"{key}: {value}" for key, value in results.items()])
        plt.text(0.1, 0.05, textstr, fontsize=12)
        plt.axis("off")
        pdf.savefig()
        plt.close()

        # If single exponential BIC is lower, plot its fit
        if bic_single < results["BIC Double Exponential"]:
            plot_results(data, lambda1, lambda2, p1, p2, lambda1_single, mean_time_single, measurement)
            pdf.savefig()
            plt.close()

        # Always plot the double exponential fit
        plot_results(data, lambda1, lambda2, p1, p2, measurement=measurement)
        pdf.savefig()
        plt.close()

        # Plot trendlines plot
        plot_results_with_trendlines(data, lambda1, lambda2, p1, p2, lambda1_single, mean_time_single, measurement)
        pdf.savefig()
        plt.close()

    print(f"Results and plots saved to {pdf_path}")


if __name__ == "__main__":
    start_time = None
    try:
        import time

        start_time = time.time()
    except:
        pass
    parser = argparse.ArgumentParser(
        prog="exponential_fit_pipeline",
        description="Perform exponential mixture model fitting on rebind-strict-boundtime.csv from the analysis pipeline."
    )
    parser.add_argument("-c", "--config", default=None, type=str, help="Path to the config TOML file")
    args = parser.parse_args()
    print(f"Running with config: {args.config}")
    main(args.config)
    if start_time is not None:
        print(f"--- {time.time() - start_time:.2f} seconds ---")
