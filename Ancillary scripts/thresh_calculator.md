# Threshold Calculator
In single-molecule tracking, particles may show apparent motion due to:
- Thermal Brownian motion (modeled by diffusion coefficient)
- Mechanical vibration/jitter from the microscope setup

To detect real molecular displacement, this module computes the step-size a particle is expected to stay within due to noise, given a desired confidence level (e.g., 99%). In another word, we want to know What is the maximum step size (step-size threshold) that can be explained by random noise (diffusion and vibration), with a given level of confidence?"
This is important because if a molecule moves more than that threshold, you can say with high confidence that it's not just noise — it's real motion.
 This module estimates the step-size threshold for Brownian motion, accounting for both diffusion and vibration-induced noise. This threshold can be used to classify molecular motion (e.g., bound vs. diffusive) in the SMol-FIESTA pipeline.
 
- Calculates expected displacement threshold using diffusion coefficient and frame interval.
- Accounts for vibration noise, modeled as a rate in pixels per second.
- Returns threshold in both pixels and micrometers.
- Supports arbitrary confidence levels (e.g., 95%, 99%).

Parameters:
- diffusion (float): Diffusion coefficient 𝐷 in µm²/s. Tells how fast a particle moves when diffusing.
- delta_t (float): Frame interval Δt in seconds.
- pixel_size (float): Pixel size in µm. Converts from µm to pixels. Needed because TrackMate and spot positions are pixel-based.
- confidence (float, optional): Confidence level (0 < confidence < 1). Default: 0.99. A value like 0.99 gives you a threshold where 99% of the expected steps (due to noise) fall below it.
- sigma_vibration (float, optional): Vibration noise rate in px/s. Default: 0.0.T. Models extra "jitter" that compounds with time.
For ecoli `1 px/s`, animal cell histone `4 px/s`, s. cerevisiae `4 px/s`

Returns:
- threshold_px (float): Step-size threshold in pixels.
- threshold_um (float): Step-size threshold in µm.
- sigma_diff (float): Diffusion-only standard deviation in pixels. 
- sigma_tot (float): Combined standard deviation in pixels (diffusion + vibration).

How it works: 
Step 1: Validate confidence value
Makes sure the confidence input is in the right range.

Step2: Compute diffusion-only standard deviation in pixels
For a particle undergoing Brownian motion with diffusion coefficient 𝐷, the standard deviation of displacement after time interval Δt is: `σ diff = 2DΔt`. Then divided by pixel_size to convert µm to pixels.

Step3:  Compute vibration noise over time
Vibration noise accumulates over time and is modeled as:  `σ vib =σ vibration × √(Δt)`. The σ vibration (vibration noise rate) of three different organism is given as an example: ecoli 1 px/s, animal cell histone 4 px/s, s. cerevisiae 4 px/s

Step4: Combine both sources of noise to get the total uncertainty (standard deviation in displacement)
Assumes the diffusion noise and vibration noise are independent, so we can add variances.

Step5: Calculate the Threshold: 
In 2D motion (like in your image plane), the step size follows a Rayleigh distribution.
Therefore, the Rayleigh distribution is used to determine the threshold value for a given confidence level (e.g. 99%): `threshold= sqrt(−2ln(1−confidence))×σ total`

For this distribution, the cumulative distribution function (CDF) is:` P(r ≤ R) = 1 - exp(-R² / (2σ²))`
We want to find R such that: `P(r ≤ R) = confidence`
Solving for R, you get: R = sqrt(-2 × σ² × ln(1 - confidence)) = c_p × σ

Therefore, c_p is a multiplier used to scale the standard deviation of noise into a step-size threshold.




