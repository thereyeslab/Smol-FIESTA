# BatchRun_Rebind.py 
This script serves as the entry point to run a series of modular analysis steps on single molecule microscopy tracking data using a configuration file in .toml format. It supports:

- Validation of file structure
- Full processing pipeline execution, including:
- Track sorting
- Cell information extraction
- Bound/diffusive classification
- Optional gap fixing
- Rebinding analysis
- Optional visualization
- Optional MSD (Mean Squared Displacement) analysis


###### Command-line Arguments:
It creates a command-line interface parser which makes it easy to run from the command line by allowing users to:
1. Specify the path to a config file (-c config.toml). If not provided, it defaults to None.
2. Run a check of file inputs (-d)
3. View a helpful message if anything is wrong or missing

It can be used like:
``` bash
python script.py -c configs/my_config.toml -d
```

###### Config file loading
It loads the config file. If no path is passed, it defaults to script-config.toml in the current directory.

###### File Check Mode
If --check-files is used:
Lists and sorts masks (.png) and CSV tracking files.
Verifies:
1. Same number of CSVs and masks.
2. That each mask has the correct number of cells (by checking the max() value in the mask).
3. That tracking data exists for each cell.

###### Main Pipeline Execution
Runs the following modules sequentially, passing the config path to each:
| Module                | Purpose                                                        |
|----------------------|----------------------------------------------------------------|
| `track_sorting`      | Reads, organizes, and filters track data (spots & tracks). |
| `cell_info`          | Analyzes mask data to extract cell-level metadata.             |
| `bound_classification` | Classifies tracks into bound or diffusive molecules using spatial features. |
| `gaps_and_fixes`     | (Optional) Fixes tracking gaps and merges trajectories.         |
| `rebind_analysis`    | Quantifies rebinding events.       |
| `visualizer`         | (Optional) Generates plots and visuals for analysis or verification. |
| `rebind_MSD`         | (Optional) Calculates MSD to characterize motion behavior.      |

Toggles for `gaps_and_fixes`, `visualizer`, and `rebind_MSD` are defined in the config file under `[toggle]`.

# track_sorting.py
This script processes TrackMate spot CSV outputs and segmented cell masks to generate cleaned, filtered, and structured track data for downstream analysis.
It expects:
- One spot CSV file per cell, per video (named with a consistent format such as `Video1_Cell_1_spotsAll.csv`)
- A corresponding mask image per video

The script performs the following steps:
1. Reads and organizes input files
2. Filters out spots outside the segmented cells
3. Groups spots by `track_id` and sorts spots within the track by frame
4. Splits tracks if there are large frame gaps
5. Filters out tracks that are too short or too long
6. Computes Euclidean distances between each spot and its neighbors in the consecutive frames (number of neighbors is configurable)
7. Optionally applies bacteria-specific cleaning for dense microscopy

The final output is a CSV file with the following columns:  
`['Video #', 'Video Name', 'Cell', 'Track', 'Frame', 'x', 'y', 'Intensity', 'Distance to Previous', 'Distance to Next']`

##### Core Functions
###### File Handling 
1. get_file_names_with_ext(path, ext) : Returns a list of file paths under a given directory with the specified extension (e.g., .csv, .png). 
2. csv_name_sort_suffix(path, suffix='spotsAll'):
Filters .csv files containing a specific suffix (default: "spotsAll") and a Cell tag in the filename.Returns a dictionary mapping each video to its corresponding cell CSV file paths.

###### parse_csv_by_mask
Reads a spot CSV file and retains only the spots that fall within a specific cell mask (based on the provided cell index).

###### track_separation(spots)
Groups spots by track_id, removes empty tracks, strips the track_id, and sorts each resulting track by frame.
Input: list of spots with [track_id, frame, x, y, quality]
Output: list of tracks, each a list of [frame, x, y, quality]

###### track_splitting_filtered 
This function splits a track into subtracks whenever there’s a gap in frames greater than gap_max, and filters out subtracks that are too short or too long (outside [len_min, len_max]). It returns: 
1. A list of filtered subtracks 
2. Number of splits made
3. Number of subtracks removed due to filtering

###### track_distance_tabulate(track, indices, dist_none)
Calculates Euclidean distances between each spot and its neighbors  (defined by indices, e.g., -1 for previous frame, +1 for next frame).
If a neighbor is not available (e.g., at the end of a track), appends dist_none.
Each spot ends up with additional values: the distances to its neighbors at the specified relative indices.

##### OPTIONAL in case the organism is small (e.g., Bacteria)
In dense microscopy images where cells are small and tightly packed, the tracking algorithm (e.g., TrackMate) can mistakenly:
- Assign a track from one cell to a neighboring cell.
- Merge two distinct tracks into a single track.

To mitigate these issues, the pipeline applies several filtering steps under the bacteria_analysis flag. These steps help detect and eliminate track artifacts caused by spatial and temporal overlaps:
1. Eliminate repeated frames within the same track:
Tracks should have one spot per frame. Multiple spots in the same frame indicate a potential mislink. Function `eliminate_repeated_frames()` keeps the most plausible spot and removes the rest.
2. Check for excessive overlap between tracks in the same frame (If merging tracks is allowed in TrackMate)
Long tracks can be cut short by the presence of random short localizations. When merging is allowed, a single track with multiple concurrent localizations will be created.
If the number of frames with multiple concurrent tracks exceeds the allowed_concurrent threshold, the entire cell's data is skipped.
3. Remove long intervals of concurrent spots
When multiple tracks overlap persistently across many frames, it may indicate sustained confusion between tracks. The `concurrent_count()` and `remove_frames()` functions detect and remove such frames if the overlap duration exceeds concurrent_max, helping to isolate reliable track segments.

##### Notes
- The number of neighbor distances to compute can be configured.
- All filtering thresholds (len_min, len_max, gap_max, allowed_concurrent, concurrent_max) are configurable.


# Cell_info.py
This script analyzes segmented cell masks (PNG images) to extract basic information about each cell:
- Number of cells in each mask
- Area of each cell (in pixels)
- Length of each cell (Euclidean diagonal of bounding box)

##### Input
- A directory of PNG mask images: {mask_path}/*.png
Each image is a labeled mask where pixel values represent individual cells (0 = background, 1 = cell 1, 2 = cell 2, etc.)
- Configuration file (.toml) specifying paths

##### Output
CSV file in the output directory with these columns: 
['Mask #', 'Mask Name', '# Cells', 'Cell', 'Area', 'Length']

# bound_classification.py 
This script performs frame-by-frame classification of particle behavior within single-molecule trajectories, using spatial proximity to neighboring spots (before and after). It assigns each spot to one of three diffusion states:
- **Freely diffusing** (0): movement is not spatially restricted.
- **Constrained diffusing** (1): movement shows intermediate confinement (relaxed vs. strict disagreement).
- **Bound** (2): movement is highly restricted based on strict distance thresholds.

**Input**: tracks.csv, which contains processed single-molecule tracks with spatial distance info between a spot and its neighbors (before and after). This is output of `track_sorting.py` module:
['Video #', 'Video Name', 'Cell', 'Track', 'Frame', 'x', 'y', 'Intensity', 'Distance to Previous', 'Distance to Next']
**Output**:
 `bound_decisions.csv`: Adds diffusion state labels to each spot. Each row corresponds to a spot and includes `original data`, along with:
- `RelaxedBound`: classification using a mild threshold
- `StrictBound`: classification using a stricter threshold
- `ConstrictedDiffusion`: 1 if relaxed and strict disagree; otherwise 0
- `Bound`: final combined score (0 = free, 1 = constrained, 2 = strictly bound)

The final columns:
[Video #, Video Name, Cell, Track, Frame, x, y, Intensity, RelaxedBound, StrictBound, ContrictedDiffusion, Bound]

##### The logic of classification:
For each spot in a track, the script evaluates its spatial proximity to neighboring spots in up to 4 previous and 4 subsequent frames. Each neighbor's distance is compared to a threshold:
- If the distance is below the threshold, it is marked as 1 (within binding range).
- If it exceeds the threshold, it is marked as 0.

The `function process_track()` generates these binary arrays for each spot under two sets of distance thresholds:
- Relaxed threshold: allows for greater distances, used for classification of constrained diffusion.
- Strict threshold: enforces tighter distance limits for binding classification.

The binary arrays are passed to two decision functions:

`determine_constrained()`: detects constrained diffusion using the relaxed threshold and requires fewer neighbors within range.
`determine_bound()`: detects strictly bound behavior using the strict threshold and requires more neighbors to be within close proximity.
Thus, the distinction between bound and constrained diffusion depends on both:
- The distance threshold used, and
- The number of neighbors within that threshold.

Finally, classification is determined as:
0 (Free): Both relaxed and strict criteria agree the spot is not bound.
2 (Bound): Both relaxed and strict criteria agree the spot is bound.
1 (Constrained diffusion): Relaxed and strict disagree (e.g., relaxed detects binding but strict does not).

# Gaps_and_fixes.py
This module further refines track classification by analyzing and correcting temporal short gaps in binding and diffusion events in single-molecule microscopy data. It also converts tracks into a format compatible with SMAUG, a software tool used to estimate diffusion coefficients.
**Main Functionalities**
1. Gap Interpolation: Gaps (missing frames between detected spots) are interpreted based on:
   - Length of the gap
   - Type of neighboring events (e.g., bound vs. diffusing)
   - Duration of surrounding events

The function `process_gaps` is responsible for that. This function creates synthetic data points for frames that were skipped during tracking (gaps), inserting them based on behavioral context using a custom rule: Gaps shorter than a configurable threshold (max_bound_gapFill) and surrounded by strict binding events can be reclassified as binding. The position of the synthetic data point is filled with [-1, -1] as a placeholder.
If a track has more than a configurable number of gaps, it will be aborted from the analysis. 
2. Event Relabeling: 
The function `event_separation`, splits a single particle track into discrete "events" segments where the particle shows consistent behavior, such as being bound, diffusive, or in transition. It also counts how often the track transitions into or out of a particular state. 
The function `pass_events`, cleans and merges behavior-labeled segments (events) based on their duration. Specifically, it:
- Relabels short events (shorter than min_time) by replacing their behavior label with a substitute.
- Merges adjacent events with the same behavior into a single event to simplify the final output.

The code performs three sequential passes of event relabeling and merging using the `pass_events` function, each targeting a different behavior type and applying specific relabeling rules based on neighboring behavior or fixed substitution.
Pass 1: Merge short free diffusion (FD) events into constrained diffusion (CD)
Pass 2: Merge SB (strict binding) into constrained diffusion (CD) based on duration
Pass 3: Merge CD into FD
The final behavior after each pass will be stored for each spot in the track. 

3. Binding Proportion Filtering (Optional)

If enabled, the script filters tracks by:
- Minimum proportion of frames spent in bound states (CD + SB)
- Maximum proportion of frames strictly bound (SB only)

Tracks outside these bounds are excluded from downstream analysis with the message **FAIL by maximum/minimum proportion binding** In the log output.  

4. Export to SMAUG Format

Tracks are translated into four SMAUG-compatible CSV files: 1. All tracks (after relabeling)
and Separated tracks by event type: 2. Diffusing, 3. Constricted (CD), and 4. Bound (SB).

##### input: 
`bound_decisions.csv` file, which is an output from the `bound_classification.py` module. 
[Video #, Video Name, Cell, Track, Frame, x, y, Intensity, RelaxedBound, StrictBound, ContrictedDiffusion, Bound]

##### Outputs
**Always generated:**
- gaps-and-fixes_decisions.csv: Final annotated track data

[Video #, Video Name, Cell, Track, Frame, x, y, Intensity, isGap, GapFixed, Pass1, Pass2, Pass3, Bound]
- SMAUG-compatible files:
    - unfiltered_spotsAll.csv
    - separated_spotsDiffusing.csv
    - separated_spotsConstricted.csv
    - separated_spotsBound.csv
    
**If filtering is active:**
- filtered_passed_spotsAll.csv ( the tracks that passed the binding proportions filtering step)
- filtered_failed_spotsAll.csv (the tracks that did not pass the binding proportions filtering step)

# rebind_analysis.py

This script identifies and analyzes rebinding events—cases where a molecule temporarily unbinds from a site (or changes its behavior) and later rebinds, either to the same location or to a different one, based on positional thresholds. It applies spatial and temporal criteria to classify a track into behavioral states and quantify their transitions over time. It works on single-molecule tracks that were previously processed for binding behavior (via bound_classification.py) and optionally gap-filled (via gaps_and_fixes.py).

This module enables:
- Detection of re-binding events, where a molecule may leave and re-enter the bound state.
- Classification of tracks into fast-diffusing, confined-diffusing, or bound behaviors.
- Quantification of behavioral transitions using transition matrices.

**Logic**: 
For each track, the function `rebind_record_proximity` detects the binding event, followed by a diffusion event, and then another binding event. Then depending on the 
- The duration of unbinding.
- The distance between the unbinding and rebinding locations.
It determines if the rebinding is valid and if the particle rebinds to the same position or another position. It also categorizes these events and tabulates statistics for further analysis.

Categorizing rebindings as:
**Same-site rebinding**: within `rebind_distance_same pixels`.
**Different-site rebinding**: beyond `rebind_distance_diff pixels`.
**Unsuccessful rebinding**: rebinding does not occur within allowed time (`max_time_rebinding`), took too long to rebind.

The script performs analysis under two levels of stringency:
- Relaxed: Allows more flexible criteria to classify an event as rebinding (e.g., "constricted diffusion" as bound).
- Strict: Requires stricter classification (e.g., only "strictly bound" state is considered bound).


--pablo--

Then, three other key functions process the behavior of the spot throughout the track. Mainely how long an event lasts and what is the next behavior and if it is an stable behavior:
`bound_record`: Records all bound state durations. 
(Bound state: Continuous stretches where the behavior matches a "bound" criterion )
What It Does:
- Tracks how long each binding event lasts.
- Filters out short, noisy events (< min_time).

Outputs:
- list of durations of all qualifying bound events

`constrained_record`:
Records constrained diffusion (where a molecule is partially confined but not fully bound) and what follows after.
(Constrained diffusion (state == 1) followed by either: Free diffusion (== 0) or Binding (== 2)
What It Does:
- Identifies segments of constrained diffusion (min_time_constrained ≤ time ≤ max_time_constrained).
- Then checks the duration of the following event . If it's long enough (≥ min_time_bound), it's considered successful.
- Outputs: Counts:
How often constrained diffusion is followed by: Diffusion (code 0) and Binding (code ≥ 2)

`diffusion_record`: Analyzes diffusive segments of a track to find how long diffusion lasts before a molecule binds again.

It records diffusion events followed by binding (based on a criteria). Evaluates if the binding is stable (duration ≥ min_time_bound) or unstable (shorter).
What It Does:
- Collects durations of each diffusion episode.
- Checks the duration of the following bound period.

Outputs:
- All diffusion times.
- A count of:
    - Short bindings (unstable)
    - Long bindings (stable)

Therefore, for each track, the script analyzes the bahavior of spot throught the track: whether there is a rebinding and if it rebinds to the same location or not, and how long is the duration of the events and what is the next event that happens. The script also calculates the proportion of different events for each track. 

Finally, for all the tracks together, the function `calculate_transition_matrices` produces a 3*3 matrix which quantifies how often transitions occur between dynamic states: fast diffusion (F.dif), confined diffusion (C.Dif), and bound (Bound). It outputs both absolute counts and transition probabilities.
**Absolute Counts**: The number of observed transitions from one state to another.
**Proportions**: The probability of transitioning from one state to another.

The final output will be saved as TrackMate output format. 

**Outputs**
All output files are saved in the directory:
``` bash
{csv_path}/{output_folder_name}/
```
| Filename | Description |
|----------|-------------|
| `RESULT_rebind.txt` | Summary of all event statistics including binding and diffusion durations, transitions, and rebinding probabilities. |
| `rebind-strict-event.csv` | Records transitions between binding events within tracks under strict criteria. |
| `rebind-strict-boundtime.csv` |Binding durations identified under strict criteria. |
| `rebind-flanked-strict-boundtime.csv` | Binding times flanked by diffusion states (i.e., D → B → D). |
| `rebind-AllDiffusion-time.csv` | all diffusion (fast + constrained) event durations. |
| `rebind-strict-rebindingtime.csv` | Time intervals between successive binding events (rebinding time). |
These files contain rebinding spot information compatible with downstream SMAUG analysis tools.



# visualizer.py
Generates RGB videos of cells with overlaid tracks from bound_decisions.csv or gaps-and-fixes_decisions.csv, highlighting spot states (e.g., Bound, Diffuse, Constricted) frame-by-frame.

**inputs**

- **Mask images** (`.png`): Segmented masks for each cell.
- **Track CSV**:
  - `bound_decisions.csv`: Raw classification output.
  - `gaps-and-fixes_decisions.csv`: Classification after interpolation and filtering.
  - **Selection controlled by** the config toggle `use_gap_fixed`.


**Logic**
- Groups tracks by video and cell.
- Initializes a full-frame background video (not cropped) with the cell body and border from the segmentation mask.
- Parses all tracks and overlays spots frame-by-frame by drawing a circular,subpixel-localized spot on the image.
- For each spot in the frame:
  - Draws a circular, subpixel-localized spot at each position.
  - Uses a distinct color per classification label (`Bound`, `Diffuse`, `Constricted`, etc.).
  - Adds a motion-blur tail between consecutive spot positions for dynamic states (e.g., Diffuse, Constricted).
- Video length is truncated to 5 frames after the last spot for conciseness.

**output**
- Mimics time-lapse videos using the tabulated spot positions and classifications.
- Videos are saved as **LZW-compressed TIFF stacks** compatible with **ImageJ**.

Videos are saved in `TYXS`: TYXS means the data is arranged with the following dimensions:
- T = Time (frames in a video)
- Y = Height (rows of pixels)
- X = Width (columns of pixels)
- S = Samples (color channels, e.g., RGB → 3 samples per pixel)
 Each frame is a 2D RGB image; background pixels are black, and molecular spots are color-coded by state.

**Notes**
- **No spatial cropping**: Ensures alignment with the original field of view and segmentation masks.
- **Tail effect** is optional and configurable (`tail_effect_enabled = true/false`).


# rebind_MSD.py
This script calculates the diffusion coefficients (D) and anomalous diffusion exponents (α) for subsections of molecular tracks, particularly focusing on segments with consistent binding state (bound/unbound). It:

- Subdivides tracks into fixed-length time windows.
- Segments these by consistent "Bound" label.
- Calculates MSD curves for each segment based on spot distances between frames.
- Fits MSD curves to extract D and α.
- Fit MSD calculations in a gaussian mixture model.
- Outputs the results (with frame info, state, and fitting values).
- Generate plots for presenting each population [bound, constrained diffusion, free diffusion] and their MSD distributions.

** Logic**
Iterates through all videos → cells → tracks and subdivides tracks into chunks of 4 frames.
Within each chunk, segments are further split based on consistent "Bound" state.
Stores each segment with: Bound state, Start/end frames, (x, y) positions. For each consistent-bound segment: Calculates MSD curve and fits using either Brownian-only or anomalous model
Returns:
D: diffusion coefficient
α: anomalous exponent (how subdiffusive or superdiffusive the motion is)

**Key Functions**
1. calculate_MSD(df)
Computes the mean squared displacement (MSD) from x/y positions. Loops over increasing time lags and averages squared displacements.

2. Fit MSD data using nonlinear least squares:
- fit_MSD: fits MSD(t) = 4D·t^α + error_correction for anomalous diffusion.
- fit_MSD_brownian: fits a simplified version with α = 1 (Brownian motion).
Both include a correction for localization error (parameter b) and camera exposure/integration timing (T_int, T_exp).

3. filter_outliers(data)
Removes statistical outliers using the IQR method. Used before log-normal fitting.

4. plot_gaussian_fit(...)
Plots histograms and log-normal fits for diffusion coefficient distributions, using log-scaled axes.

**Input**:
    Either of:
- gaps-and-fixes_decisions.csv
- bound_decisions.csv
Determined by parameter: {use_gap_fixed}

**Ouput**
(Check the input/output documentation for more detail.)
`Diffusion_Coefficient_Calculation.csv` : each row in the output CSV file contains:
- Segment identity (Video, Cell, Track, Event)
- MSD-based measurements:
- Diffusion Coefficient (D)
- Alpha (diffusion exponent)
- Metadata:
    - Bound state
    - Initial Frame, Last Frame
    - Frame counts

`Diffusion_Coefficient_Plots.pdf`:
Plot of distribution of apparent diffusion coefficient for all behaviours with Gaussian Fits.




