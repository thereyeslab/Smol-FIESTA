#!/usr/bin/env python3
"""
threshold_calculator.py

Module for calculating diffusion‐based step‐size thresholds with vibration noise
that compounds with the frame interval.

Usage (from another script):
    from threshold_calculator import calculate_threshold

    # diffusion coefficient D (µm²/s), frame interval Δt (s), pixel size (µm)
    D = 0.012
    dt = 0.01
    px = 0.130

    # vibration noise rate σ_vib = 5 px/s
    thresh_px, thresh_um, sig_diff, sig_tot = calculate_threshold(
        D, dt, px, confidence=0.99, sigma_vibration=5.0
    )
"""

import math

def calculate_threshold(
    diffusion: float,
    delta_t: float,
    pixel_size: float,
    confidence: float = 0.99,
    sigma_vibration: float = 0.0
) -> tuple[float, float, float, float]:
    """
    Calculate the step‐size threshold for Brownian diffusion with vibration noise
    that compounds with the interval.

    Parameters
    ----------
    diffusion : float
        Diffusion coefficient D in µm²/s.
    delta_t : float
        Frame interval Δt in seconds.
    pixel_size : float
        Pixel size in µm.
    confidence : float, optional
        Desired confidence level (0 < confidence < 1). Default is 0.99.
    sigma_vibration : float, optional
        Vibration noise **rate** in pixels per second (px/s). Default 0.0.

    Returns
    -------
    threshold_px : float
        Threshold step‐size in pixels.
    threshold_um : float
        Threshold step‐size in µm.
    sigma_diff : float
        Diffusion‐only σ in pixels.
    sigma_tot : float
        Combined σ (diffusion + vibration) in pixels.
    """
    if not (0.0 < confidence < 1.0):
        raise ValueError("confidence must be between 0 and 1 (exclusive)")

    # Diffusion σ (µm) → px
    sigma_diff = math.sqrt(2 * diffusion * delta_t) / pixel_size

    # Vibration σ accumulates as σ_vib_rate * sqrt(Δt)
    sigma_vib_px = sigma_vibration * math.sqrt(delta_t)

    # Total noise
    sigma_tot = math.sqrt(sigma_diff**2 + sigma_vib_px**2)

    # Rayleigh percentile multiplier
    c_p = math.sqrt(-2.0 * math.log(1.0 - confidence))

    threshold_px = c_p * sigma_tot
    threshold_um = threshold_px * pixel_size

    return threshold_px, threshold_um, sigma_diff, sigma_tot


if __name__ == "__main__":
    # Example parameters
    D = 0.008         # µm²/s diffusion coeffiecient of bound
    dt = 0.01         # s (10 ms)
    px_size = 0.094   # µm
    sigma_vib = 4 # px/s (vibration noise rate) ecoli 1, animal cell histone 4, s. cerevisiae 4
    conf = 0.99       # 99% confidence

    thresh_px, thresh_um, sig_diff, sig_tot = calculate_threshold(
        D, dt, px_size, confidence=conf, sigma_vibration=sigma_vib
    )

    print(f"σ_diff       = {sig_diff:.3f} px")
    print(f"σ_vibration  = {sigma_vib:.3f} px/s → {sigma_vib*math.sqrt(dt):.3f} px over Δt")
    print(f"σ_total      = {sig_tot:.3f} px")
    print(f"threshold    = {thresh_px:.3f} px ({thresh_um:.3f} µm) at {conf*100:.1f}%")
