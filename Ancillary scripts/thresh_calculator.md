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
**Step 1**: Validate confidence value
Makes sure the confidence input is in the right range.

**Step2**: Compute diffusion-only standard deviation in pixels
For a particle undergoing Brownian motion with diffusion coefficient 𝐷, the standard deviation of displacement after time interval Δt is: `σ diff = 2DΔt`. Then divided by pixel_size to convert µm to pixels.

**Step3**:  Compute vibration noise over time
Vibration noise accumulates over time and is modeled as:  `σ vib =σ vibration × √(Δt)`. The σ vibration (vibration noise rate) of three different organism is given as an example: ecoli 1 px/s, animal cell histone 4 px/s, s. cerevisiae 4 px/s

**Step4**: Combine both sources of noise to get the total uncertainty (standard deviation in displacement)
Assumes the diffusion noise and vibration noise are independent, so we can add variances.

**Step5**: Calculate the Threshold: 
In 2D motion (like in your image plane), the step size follows a Rayleigh distribution.
Therefore, the Rayleigh distribution is used to determine the threshold value for a given confidence level (e.g. 99%): `threshold= sqrt(−2ln(1−confidence))×σ total`

For this distribution, the cumulative distribution function (CDF) is:` P(r ≤ R) = 1 - exp(-R² / (2σ²))`
We want to find R such that: `P(r ≤ R) = confidence`
Solving for R, you get: R = sqrt(-2 × σ² × ln(1 - confidence)) = c_p × σ

Therefore, c_p is a multiplier used to scale the standard deviation of noise into a step-size threshold.

#### How to Use the Script

**1. As a standalone test:**
Open a terminal and run:
```bash
python threshold_calculator.py
```
This runs the default example and prints the threshold in pixels and microns. You can change the values of the default example as desired and run the script to get the threshold proper for your own data. 
You need to change this part of the code: 
```bash
D = 0.008         # µm²/s diffusion coeffiecient of bound
dt = 0.01         # s (10 ms)
px_size = 0.094   # µm
sigma_vib = 4 # px/s (vibration noise rate) ecoli 1, animal cell histone 4, s. cerevisiae 4
conf = 0.99  
````

**2. From another script:**
```bash
from threshold_calculator import calculate_threshold

# Define your experiment parameters:
D = 0.012       # diffusion coefficient (µm²/s)
dt = 0.01       # frame interval (s)
px = 0.130      # pixel size (µm)
vibration = 5.0 # vibration noise rate (px/s)

# Get threshold and noise estimates
threshold_px, threshold_um, sigma_diff, sigma_total = calculate_threshold(
    diffusion=D,
    delta_t=dt,
    pixel_size=px,
    confidence=0.99,
    sigma_vibration=vibration
)

print("Step-size threshold:", threshold_px, "px (", threshold_um, "µm )")
```



#### Troubleshooting

| **Issue**                                                    | **Possible Cause**                                                                                  | **Solution**                                                                                                         |
| ------------------------------------------------------------ | --------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `ValueError: confidence must be between 0 and 1 (exclusive)` | Confidence value passed is ≤ 0 or ≥ 1                                                               | Ensure `confidence` is a float like `0.95` or `0.99`, not `1.0` or `100`                                             |
| `ZeroDivisionError` or `inf threshold_px`                    | Pixel size (`px`) is set to 0                                                                       | Make sure `pixel_size` is a positive number in µm/pixel                                                              |
| Output threshold seems unrealistically high or low           | Incorrect units or values for `D`, `dt`, or `px`                                                    | Double-check that:<br>- `D` is in µm²/s<br>- `dt` is in seconds<br>- `px` is in µm                                   |
| `NameError` for `calculate_threshold`                        | You didn’t import the function in your script                                                       | Add `from threshold_calculator import calculate_threshold`                                                           |
| Script runs but prints nothing                               | You're running the module in an interactive notebook or haven’t defined `if __name__ == "__main__"` | To teststandalone: run directly as `python threshold_calculator.py` or copy the example usage into a script         |
| Confused about how to apply threshold                        | The script only calculates the cutoff; it doesn't apply it to real data                             | Use the `threshold_px` output as a cutoff for step sizes in your trajectory data (in the `Bound-Classification` section in the config file, you can set the value of  `distance_threshold_strict` using the given threshold.    |
| Unsure what “σ\_vibration” means                             | You don’t have a good estimate for vibration noise                                                  | Use empirical estimates:<br>- E. coli: `σ_vib ≈ 1 px/s`<br>- Yeast or histone in mammalian cells: `σ_vib ≈ 4–5 px/s`  as a stat. Adjust if needed.|

