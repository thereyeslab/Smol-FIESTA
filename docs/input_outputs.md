# batchRun_rebind.py
Master script that runs the full single-molecule analysis pipeline using a specified TOML config file. It orchestrates multiple modules for rebinding analysis. 

**Inputs**
1. Config File (.toml)
Specifies paths, parameters, and module toggles ( check the `config.md` for details)

2. TrackMate spots files
Format: one CSV per cell per video.
Each row corresponds to a detected spot. Must follow this naming convention: `*_Cell_*_spotsAll.csv` (ex: `VideoName_Cell_1_spotsAll.csv`)
Typical columns include:
    | Column Name   | Description                              |
    |---------------|------------------------------------------|
    | TRACK_ID      | ID of the track the spot belongs to.     |
    | POSITION_T    | Time point (frame number).               |
    | POSITION_X    | X coordinate of the spot (in pixels).    |
    | POSITION_Y    | Y coordinate of the spot (in pixels).    |
    | INTENSITY     | Intensity value of the spot.             |

3. Segmentation Masks (.png)
Grayscale images where each cell has a unique integer label (0 = background). One mask per video.

**Note** : The path to the spot files and mask files should be stored in the config file. 

**Outputs**: it will be discussed for each step of the pipeline bellow. 

# Track_sorting.py
This script processes TrackMate spot CSV files and segmentation masks and outputs the cleaned and structured track data for downstream analysis.(see documentation.md for more details)

**Inputs**:
Since it is the first step of the pipline, the inputs are as in the master script (`batchRun_rebind.py`).
1. Config File (.toml)
Specifies paths, parameters, and module toggles ( check the `config.md` for details)
2. TrackMate spots files
Format: one CSV per cell per video.
Each row corresponds to a detected spot. Must follow this naming convention: `*_Cell_*_spots.csv` (ex: `VideoName_Cell_1_spots.csv`)
Typical columns include:
    | Column Name   | Description                              |
    |---------------|------------------------------------------|
    | TRACK_ID      | ID of the track the spot belongs to.     |
    | POSITION_T    | Time point (frame number).               |
    | POSITION_X    | X coordinate of the spot (in pixels).    |
    | POSITION_Y    | Y coordinate of the spot (in pixels).    |
    | INTENSITY     | Intensity value of the spot.             |

3. Segmentation Masks (.png)
Grayscale images where each cell has a unique integer label (0 = background). One mask per video.

**Parameters** (from script-config.toml)

Under `[track-sorting]`:
| Parameter                  | Description                                                                 |
| -------------------------- | --------------------------------------------------------------------------- |
| `allowed_gap_max`          | Max frame gap allowed within a track before it’s split.                     |
| `allowed_track_length_min` | Minimum allowed track length.                                               |
| `allowed_track_length_max` | Maximum allowed track length.                                               |
| `dist_range`               | Number of frames before and after the spot to use for distance calculation.          |
| `allowed_overlap`          | allowed total repeated frames in tracks in a cell, eliminate all tracks if exceeds                        |
| `concurrent_max`           | Max number of concurrent frames in *different* tracks in a cell, eliminate the frames if exceeds    |

Under `[toggle]`
| Parameter           | Description                                                      |
| ------------------- | ---------------------------------------------------------------- |
| `bacteria_analysis` | If True, additional filters are applied for small-cell analysis. |

**Outputs**
1. tracks.csv:
File: {csv_path}/{output_folder_name}/Intermediates/tracks.csv
A cleaned and annotated DataFrame where each row represents a spot in a valid track, with columns:
Video #, Video Name, Cell, Track, Frame, x, y, Intensity

Dist -{n} to Dist +{n}: Spot distances to other frames within the specified dist_range.

2. script-config.toml:
A copy of the configuration file, saved alongside the output CSV for reproducibility.

# cell_info.py
This script computes basic information about each cell in segmentation mask images.
**Inputs**
Segmentation Masks (.png)
Grayscale images where each cell has a unique integer label (0 = background). One mask per video.

**Outputs**
_cell-info.csv
| Column Name   | Description                                            |
| ------------- | ------------------------------------------------------ |
| `Mask #`    | Index of the mask in the input list                    |
| `Mask Name` | File name of the mask image                            |
| `# Cells`   | Total number of segmented cells in the mask            |
| `Cell`      | Index of the cell within the mask (starting from 1)    |
| `Area`      | Pixel count (area) of the individual cell              |
| `Length`    | Euclidean diagonal of the bounding box around the cell |

Output format: {csv_path}/{output_folder_name}/_cell-info.csv

# bound_classification.py
This script analyzes parsed tracks to determine the binding behavior of particles frame-by-frame using distance thresholds. It produces a labeled decision table indicating whether each particle is freely diffusing, constrained, or bound.

**Inputs**
tracks.csv that is outputed from the `track_sorting.py`
{csv_path}/{output_folder_name}/Intermediates/tracks.csv

**Parameters**
| **Parameter**               | **Description**                                                                                                     |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `distance_threshold`        | Distance in pixels within which two spots must fall to be considered **constrained** (used for `RelaxedBound = 1`). |
| `distance_threshold_strict` | Stricter distance threshold (in pixels) for two spots to be considered **strictly bound** (`StrictBound = 1`).      |

**Outputs**
format: {csv_path}/{output_folder_name}/Intermediates/bound_decisions.csv

bound_decisions.csv
Original spot data with the following columns added: 
| **Column Name**        | **Type** | **Description**                                                                          |
| ---------------------- | -------- | ---------------------------------------------------------------------------------------- |
| `RelaxedBound`         | Integer  | 1 if constrained diffusion under the relaxed threshold, else 0                           |
| `StrictBound`          | Integer  | 1 if bound under the strict threshold, else 0                                            |
| `ConstrictedDiffusion` | Integer  | 1 if constrained but not strictly bound (i.e., `RelaxedBound` XOR `StrictBound`), else 0 |
| `Bound`                | Integer  | Final state classification: 0 = free, 1 = constrained, 2 = bound                 |

# Gaps_and_fixes.py
This script further refines track behavior interpretation by
- interpolating the gap behavior based on neighboring events and their durations.
- introducing conditional filtering by overall binding event probability.


**Input**: 
bound_decisions.csv outputted from the bound_classification.py
Format: {csv_path}/{output_folder_name}/Intermediates/bound_decisions.csv

**Parameters**
| **Parameter**          | **Type** | **Description**                                                                       |
| ---------------------- | -------- | ------------------------------------------------------------------------------------- |
| `min_time_strict`      | Integer  | Minimum number of frames required to classify an event as strictly bound.             |
| `min_time_constrained` | Integer  | Minimum number of frames for an event to be considered constrained diffusion.         |
| `min_time_diffusion`   | Integer  | Minimum number of frames to classify an event as free diffusion.                      |
| `max_bound_gapFill`    | Integer  | Maximum gap (in frames) that can be interpolated as part of a binding event.          |
| `min_prop_binding`     | Float    |  Minimum proportion of frames spent in bound states (used for filtering).      |
| `max_prop_binding`     | Float    | Maximum proportion of frames strictly bound (used to remove overly bound tracks). |

Conditional Parameters (used for advanced filtering)
| **Parameter**        | **Type** | **Description**                                                                     |
| -------------------- | -------- | ----------------------------------------------------------------------------------- |
| `allowed_gapTotal`   | Integer  | Maximum number of total gap frames allowed to be interpolated in a single track.    |
| `Prop_check_over`    |     |        |
| `Ratio_turnover_dtf` |     |  |
| `Ratio_turnover_all` |     |                                   |

**Ouputs**
 `gaps-and-fixes_decisions.csv` in the {csv_path}/{output_folder_name}/Intermediates/gaps-and-fixes_decisions.csv
Contains the original spot data with the added columns below:
| **Column Name** | **Description**                                                                                    |
| --------------- | -------------------------------------------------------------------------------------------------- |
| `isGap`         | 1 if this frame is a gap (interpolated, not an observed spot); otherwise 0.                        |
| `GapFixed`      | 1 if this gap was filled during the interpolation process; otherwise 0.                            |
| `Pass1`         | Behavior label after first event relabeling pass: Fix short free diffusion (FD) → CD.                   |
| `Pass2`         |Behavior label after second event relabeling pass: Fix short SB → CD.                |
| `Pass3`         | Behavior label after third event relabeling pass: Fix short CD → FD.         |
| `Bound`         | Final state assigned to the spot|

**Note**
The Pass1, Pass2, and Pass3 columns store the intermediate behavior state at each relabeling pass. This provides traceability of how a spot's label evolved across steps.


# rebind_analysis.py
This script identifies and characterizes rebinding events in single-molecule tracks. It analyzes the dynamics of particles transitioning from bound states to free states and back to binding. It also evaluates whether the rebindings occur to the same molecular site or to different ones, based on spatial displacement.

**input**
One of the following files (based on the use_gap_fixed flag):
        {csv_path}/{output_folder_name}/Intermediates/gaps-and-fixes_decisions.csv --> outputted from `Gaps_and_fixes.py`
        {csv_path}/{output_folder_name}/Intermediates/bound_decisions.csv --> outputted from `bound_classification.py`
        
**Parameters**
| Parameter                       | Description                                                                              |
| ------------------------------- | ---------------------------------------------------------------------------------------- |
| `rebind_distance_same`          | Max spatial distance (pixels) to classify a rebind as returning to the **same** site.    |
| `rebind_distance_diff`          | Min distance to classify as rebinding to a **different** site.                           |
| `min_time_bound_strict`         | Minimum duration (frames) for a valid **strict binding** event.                          |
| `min_time_bound_constrained`    | Minimum duration (frames) for a valid **constrained diffusion** event.                   |
| `min_time_rebinding_relaxed`    | Minimum duration (frames) for rebinding under **relaxed** rules.                         |
| `min_time_rebinding_strict`     | Minimum duration (frames) for rebinding under **strict** rules.                          |
| `min_time_diffusion`            | Minimum duration (frames) required for the diffusion segment between two binding events. |
| `min_time_diffusion_subsequent` | Minimum time required after diffusion to qualify as a rebinding event.                   |
| `max_time_rebinding`            | Maximum allowed time (frames) between unbinding and rebinding to count as valid.         |
| `max_time_constrained`          | Similar to above, but for constrained state; exceeding this cancels rebinding status.   |


**ouputs**
1. `RESULT_rebind.csv` : Summary of rebinding statistics (counts, proportions).


**Rebinding section**
| Metric                              | Description                                      |
| ----------------------------------- | ------------------------------------------------ |
| Rebinding (Events)                  | Number of successful vs. failed rebinding events + (Probability of Success) Ratio of successful rebindings to total attempts|
**Note**
Unsuccessful rebinding: rebinding does not occur within allowed time (`max_time_rebinding`), took too long to rebind.

**All Events section**
| Metric                        | Description                            |
| ----------------------------- | -------------------------------------- |
| Fast Diffusion       | number and proportion of events in fast diffusive state  |
| Constrained Diffusion| number and proportion of events in constrained diffusive state    |
| Bound                | number and proportion of events in bound state |

2. `RESULT_rebind.txt`
| Section                     | Column Name                        | Description                                                                          |
| --------------------------- | ---------------------------------- | ------------------------------------------------------------------------------------ |
| **Event Mean Time**         | `Bound`                            | Duration the molecule remained in the bound (immobile) state             |
|                             | `C. Diffusion`                     | Duration in the confined diffusion state (restricted motion)             |
|                             | `F. Diffusion`                     | Duration in the fast diffusion (free movement) state                     |
|                             | `Diffusion Time`                   | Combined time in any diffusion state (Confined + Fast)                   |
|                             | `Rebind Time`                      | Time (in frames) it took for a molecule to return to the bound state after unbinding |
|                             | `count`                            | Number of trajectories or events contributing to the statistics                      |
|                             | `mean`, `std`, `min`               | Mean, standard deviation, and minimum duration (frames) for that event type          |
|                             | `25%`, `50%`, `75%`                | First, second (median), and third quartile durations                                 |
|                             | `max`                              | Maximum observed duration in frames                                                  |
| **Rebinding Probability**   | `Successful`                       | Number of unbound molecules that re-entered the bound state within 50 frames         |
|                             | `Unsuccessful`                     | Number of unbound molecules that did **not** rebind within 50 frames                 |
|                             | `Probability`                      | Rebinding success rate = Successful / (Successful + Unsuccessful)                    |
| **Event Proportions**       | `Event`                            | The type of motion state: Bound, Confined Diffusion, or Fast Diffusion               |
|                             | `Count`                            | Number of sub-tracks or events assigned to each motion state                         |
|                             | `Proportion`                       | Fraction of total events that belong to each state                                   |
| **Transition Matrix (3x3)** | `From_X.To_Y`                      | Probability per frame of switching from one state (X) to another state (Y)           |
|                             | `From_F.dif`                       | Transitions starting from Fast Diffusion                                             |
|                             | `From_C.Dif`                       | Transitions starting from Confined Diffusion                                         |
|                             | `From_Bound`                       | Transitions starting from Bound                                                      |
|                             | `To_F.dif`, `To_C.Dif`, `To_Bound` | Probabilities of entering each state from the current one                            |
| **Transition Matrix (2x2)** | `From Diffusion`                   | Combined transitions from both diffusion states (Fast + Confined) to Bound or itself |
|                             | `From Bound`                       | Transitions from Bound to Diffusion or itself                                        |
|                             | `To Diffusion`                     | Probability of transitioning into any diffusion state                                |
|                             | `To Bound`                         | Probability of transitioning into Bound state                                        |


3. `rebind-Events.csv`
This file records individual state transitions per track, capturing detailed rebinding events.

| Column Name | Description                                                                         |
| ----------- | ----------------------------------------------------------------------------------- |
| `Video #`   | Index or ID of the video in the dataset                                             |
| `Cell`      | Identifier for the cell in the video (if multiple cells are tracked)                |
| `Track`     | Unique identifier for a single molecule’s trajectory                                |
| `Event`     | Index of the event (e.g., first event = 0, second event = 1, etc.) within the track       |
| `time`      | Duration of this event — how long the molecule stayed in this state |
| `type`      | Type of event/state (see below)                                                     |

type values (event categories):

| Type            | Description                                                               |
| --------------- | ------------------------------------------------------------------------- |
| `Bound`         | Molecule remained in the immobile (bound) state                           |
| `C.Diffusion`   | Confined diffusion: molecule shows restricted movement                    |
| `FastDiffusion` | Free/fast diffusion: molecule moves relatively freely                     |
| `SearchTime`    | Time between unbinding and next rebinding (i.e., unbound search interval) |
| `Rebinding`     | Event marking the return to the bound state after unbinding               |



# Visualizer.py
Runnable Script if run as __main__
Generate videos with labelled tracks from bound-decisions or gap-and-fixes.
- Mimic timelapse videos from the tabulated spot positions and their behavior.
- Spots for each behavior are drawn using a single color per class.
- A motion-blur tail is drawn for moving molecules (classes 0 and 1) to indicate motion.
- Video is truncated to 5 frames after the last spot.
- Output file is named using the mask name.

**input**
One of the following files (based on the use_gap_fixed flag):
        {csv_path}/{output_folder_name}/gaps-and-fixes_decisions.csv --> outputted from `Gaps_and_fixes.py`
        {csv_path}/{output_folder_name}/bound_decisions.csv --> outputted from `bound_classification.py`
**Parameters**: 
`allowed_track_length_max`: Max number of frames in the video (truncates longer sequences).

`use_gap_fixed` (boolean):
Whether to use output from gaps_and_fixes.py (True) or bound_classification.py (False).

**Visualizer-specific** (under [visualizer] in config):
`spot_diameter`: Diameter (in pixels) of drawn molecular spots.

`tail_effect_enabled`: If True, diffusing molecules (classes 0 and 1) will have a motion blur tail indicating their path.



**Output**
{csv_path}/{output_folder_name}/Visualizer/{mask_name}.tif
| Class/Behavior | Value in CSV | Color   | Drawn Tail? |
| -------------- | ------------ | ------- | ----------- |
| Diffusive        | `0`          | Magenta | ✅           |
| Constricted    | `1`          | Blue    | ✅           |
| Bound          | `2`         | Green   | ❌           |
| Gap/Missing    | NaN coords   | Black   | ❌           |

# rebind_MSD.py

Performs MSD calculations with gaussian mixture model.
- Calculate MSD based on spot distances between frames.
- Fit MSD calculations in a gaussian mixture model.
- Generate plots for presenting each population [bound, constrained diffusion, free diffusion] and their MSD distributions.

**input**
One of the following files (based on the use_gap_fixed flag):
        {csv_path}/{output_folder_name}/gaps-and-fixes_decisions.csv --> outputted from `Gaps_and_fixes.py`
        {csv_path}/{output_folder_name}/bound_decisions.csv --> outputted from `bound_classification.py`

** parameters**

**From Config**

| Parameter           | Description                                         |
|---------------------|-----------------------------------------------------|
| `pixel_size_um`     | Converts pixel coordinates to microns               |
| `time_interval_ms`  | Time between frames in milliseconds                                |
| `use_brownian_only` | If true, force α = 1 (pure Brownian motion)         |

**in Script**

| Parameter              | Value   | Description                                      |
|------------------------|---------|--------------------------------------------------|
| `frames_per_subdivision` | 4       | MSD calculated in 4-frame segments               |
| `b`                    | 0.0     | Localization error in microns (for correcting MSD) |
| `T_int`                | 0.01    | Interval time (seconds)                       |
| `T_exp`                | 0.01    | Exposure time (seconds)                          |

**Diffusion Coefficient Filtering**

| Parameter            | Value   | Description                     |
|----------------------|---------|---------------------------------|
| `cutoff_bound`       | 0.02    | Max bound diffusion coefficient |
| `cutoff_cdiffusion`  | 2.5     | Max confined diffusion  coefficient        |
| `cutoff_fdiffusion`  | 5       | Max free diffusion  coefficient            |
| `min_bound`          | 0.0001  | Min bound diffusion coefficient |
| `min_cdiffusion`     | 0.001   | Min confined diffusion          |
| `min_fdiffusion`     | 0.001   | Min free diffusion              |

**Outputs**
1. `Diffusion_Coefficient_Calculation.csv` in {csv_path}/{output_folder_name}/Intermediates/Diffusion_Coefficient_Calculation.csv

| Column Name              | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| `Video #`                | Index or identifier of the video being analyzed                             |
| `Video Name`             | Filename or name of the video                                               |
| `Cell`                   | Identifier for the cell in the video (if multiple cells are tracked)        |
| `Track`                  | Unique track ID for a particle/trajectory                                   |
| `Event`                  | Type of event or track annotation (e.g., entry/exit event, interaction)     |
| `Diffusion Coefficient` | Estimated diffusion coefficient (µm²/s) from MSD fit                         |
| `Alpha`                  | Anomalous diffusion exponent (α); indicates the diffusion regime             |
| `Raw_Count`              | Total number of frames in the sub-track                                    |
| `Valid_Count`            | Number of valid data points used in the MSD fitting                         |
| `Bound`                  | Classification of the particle (e.g., Bound, CD, Unbound) based on D/α    |
| `Initial Frame`          | Index of the first frame of the sub-track                                |
| `Last Frame`             | Index of the last frame of the sub-track                                     |

2. PDF Plot:
{csv_path}/{output_folder_name}/Diffusion_Coefficient_Plots.pdf
Contains:
    - Individual histograms for each population: Bound (green), Constrained (blue), Free (red)
    - Combined histogram with overlaid Gaussian fits (optional)

# TwoExponential.py

This module fits dwell-time data (bound and search events) to exponential models using either single- or double-exponential fitting. It classifies events using posterior probabilities and calculates kinetic transition probabilities between states.

**Inputs**
CSV file: rebind-Events.csv
(Generated by bound_classification.py or gaps_and_fixes.py)
| Column    | Description                                                      |
| --------- | ---------------------------------------------------------------- |
| `Video #` | ID or index of the video                                         |
| `Cell`    | Cell identifier (if multiple cells per video)                    |
| `Track`   | Unique trajectory ID for a single molecule                       |
| `Event`   | Index of the event within the trajectory (0 = first event, etc.) |
| `time`    | Duration (in frames) that the molecule stayed in a state         |
| `type`    | Type of event (e.g., `Bound`, `SearchTime`)                      |

**Outputs**
1. `Exp1_probabilities.csv`
Posterior probabilities for each dwell event, classified under the first exponential component.

| Column Name           | Description                                                   |
| --------------------- | ------------------------------------------------------------- |
| `Video #`             | Identifier of the input video                                 |
| `Cell`                | Cell ID if applicable                                         |
| `Track`               | Track ID of the molecule                                      |
| `Event`               | Index of the dwell event                                      |
| `StartFrame`          | Frame where the event begins                                  |
| `EndFrame`            | Frame where the event ends                                    |
| `time`                | Duration of the event in seconds                              |
| `type`                | Event type (`Bound` or `SearchTime`)                          |
| `Posterior_Prob_Exp1` | Posterior probability that the event belongs to Exponential 1 |

2. `transition_probabilities.csv`
Transition matrix estimating the probability of moving between kinetic states (e.g., Bound_Exp2 → SearchTime_Exp1). The transition matrix is a square matrix where each entry (i, j) gives the probability of transitioning from state i to state j.

3. `2Exponential_test.pdf`
The PDF output includes fit parameters and model diagnostics for bound and search times. A summary figure that includes histogram plots of dwell time or search time durations and their exponential fits, along with statistical annotations.Below are the descriptions of each reported parameter:

| Field                                      | Description                                                                                                                                                                                                              |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Estimated Mean Time 1 (non-corrected)**  | The mean dwell time (in seconds) associated with the first exponential component, typically representing short binding events. Calculated as the inverse of the decay rate (1/λ₁).                                       |
| **Estimated Decay Rate 1 (non-corrected)** | The decay rate (λ₁) for the first exponential component. Higher λ implies faster decay (shorter events).                                                                                                                 |
| **Estimated Mean Time 2 (non-corrected)**  | The mean dwell time associated with the second exponential component, typically representing longer binding events.                                                                                                      |
| **Estimated Decay Rate 2 (non-corrected)** | The decay rate (λ₂) for the second exponential component. Lower λ implies slower decay (longer events).                                                                                                                  |
| **Fraction of Exp Component 1**            | The proportion of all events that are assigned to the first exponential component (short-lived binding).                                                                                                                 |
| **Fraction of Exp Component 2**            | The proportion of events assigned to the second exponential component (long-lived binding).                                                                                                                              |
| **Weight Fraction 1**                      | The contribution of the first exponential component to the total fit, weighted by both its fraction and decay characteristics.                                                                                           |
| **Weight Fraction 2**                      | The contribution of the second exponential component to the total fit.                                                                                                                                                   |
| **BIC Double Exponential**                 | The Bayesian Information Criterion for the double-exponential model. A lower BIC indicates a better fit while penalizing for added complexity (i.e., more parameters). Useful for comparing to a single-exponential fit. |

