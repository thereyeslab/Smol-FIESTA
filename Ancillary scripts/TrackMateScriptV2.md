# TrackMateScriptV2.py

This script automates tracking of single molecules or particles in a set of .tif microscopy images using TrackMate (Fiji plugin). IIt is designed to batch-process .tif microscopy videos, optionally with segmentation outlines, and extract per-track and per-spot features, saving them as .csv files for downstream analysis.

This is a standard script used for tracking and also could be used to output the track and spot files with the patern compatible to use in the SMOL pipeline. 

##### High-Level Workflow
The script performs the following main tasks:
- Loads configuration parameters from a .txt config file (Parse tracking, filtering, ROI, and file path settings from a config file).

- Loads .tif image files for analysis and their corresponding ROI outlines (if segmentation is enabled).

- For each image, perform segmentation (if enabled). Initializes tracking parameters for each video, runs tracking using TrackMate with specified detection and linking settings.

- Optionally segments ROIs (cells). For each ROI/cell/video, runs TrackMate with defined parameters: 1.Detects spots (via DoG, LoG, or thresholding) 2. Tracks molecules using the [LAP tracker](https://imagej.net/plugins/trackmate/trackers/lap-trackers). 3. Applies user-defined filters.

- Extracts features (speed, intensity, location, etc.) from each track and spot.

- Saves results to per-cell CSV files (*_tracks.csv, *_spots.csv) in a structured output directory.

#### Input Files
1. TIFF Videos (Must be in .tif format)
Located in a directory defined in the config file under [Paths] -> dir_tiff.

2. ROI Outline Files (Optional)
Used if segmentation_enabled = True.
Each .txt file contains 1+ polygonal ROI per frame. It is potentially created by [Cellpose2](https://cellpose.readthedocs.io/en/latest/).Should be in the same order as the TIFFs and contain comma-separated X,Y coordinates for each ROI.

3. Config File
A text file like Trackmate_configFile.txt containing required parameters.

#### How It Works
Step-by-Step Pipeline:
1. main() is called:
    - Loads the config.
    - Scans the TIFF and outline directories.
    - Verifies file matches (if segmentation is enabled).
    - Prepares an output folder: AnalysisSMolFiesta.

2. For each TIFF file:
If segmentation is on, iterates through all ROIs (cells) in the corresponding .txt file.
For each ROI or whole image:
    - Configures TrackMate detector, tracker, and filters.
    - Runs tracking with a timeout safety check (60s per ROI).
    - Extracts features: speed, intensity, coordinates, etc.
    - Saves results as:
        - *_Cell_<id>_tracks<suffix>.csv
        - *_Cell_<id>_spots<suffix>.csv


#### Key Modules & Functions
| Function             | Purpose                                                     |
| -------------------- | ----------------------------------------------------------- |
| `main()`             | Orchestrates the pipeline, loads files and config           |
| `load_config(path)`  | Parses the INI-style config file                            |
| `tracking()`         | Runs tracking per video and ROI                             |
| `configure()`        | Sets up TrackMate detector, tracker, filters                |
| `run_with_timeout()` | Ensures no infinite loops during tracking                   |
| `track_process()`    | Runs TrackMate, extracts features, writes output            |
| `load_rois()`        | Loads polygon ROIs from text files, optionally expands them |


#### Outputs
Saved to: <dir_tiff>/AnalysisSMolFiesta/

For each ROI (or full image):
- Tracks CSV: one row per track with summary stats: Mean speed, quality, duration, average location, etc:

| Column Name | Description                                          |
| ----------- | ---------------------------------------------------- |
| `Track_IDs` | Unique ID for each detected track.                   |
| `spt_tr`    | Number of spots in the track.                        |
| `spt_widt`  | Average radius/width of the spots in the track.      |
| `mean_sp`   | Mean speed of spots within the track.                |
| `max_sp`    | Maximum speed within the track.                      |
| `min_sp`    | Minimum speed within the track.                      |
| `med_sp`    | Median speed within the track.                       |
| `std_sp`    | Standard deviation of speed in the track.            |
| `mean_q`    | Mean quality score of the track (TrackMate feature). |
| `max_q_tr`  | Maximum spot quality within the track.               |
| `min_q_tr`  | Minimum spot quality within the track.               |
| `med_q_tr`  | Median spot quality within the track.                |
| `std_q_tr`  | Standard deviation of spot quality within the track. |
| `tr_dur`    | Duration of the track in frames.                     |
| `tr_start`  | Start frame of the track.                            |
| `tr_fin`    | End frame of the track.                              |
| `x_lc`      | Average X-coordinate of the track’s center of mass.  |
| `y_lc`      | Average Y-coordinate of the track’s center of mass.  |

- Spots CSV: one row per spot: Track ID, frame number, x/y location, intensity, radius.

| Column Name   | Description                                                |
| ------------- | ---------------------------------------------------------- |
| `tr_identifi` | Track ID to which the spot belongs.                        |
| `tr_fram`     | Frame number where the spot was detected.                  |
| `x_tr`        | X-coordinate of the spot.                                  |
| `y_tr`        | Y-coordinate of the spot.                                  |
| `inten`       | Total intensity of the spot. |



### Dependencies
- Python 3.x

- configparser, csv, os, shutil, math, time, threading

- ImageJ/Fiji Java integration (e.g., via PyImageJ or running inside Jython environment in Fiji)

**IMPORTANT NOTE**: This script is meant to be run inside Fiji or using a Fiji-compatible Python interpreter that supports Java-based modules like TrackMate, ImagePlus, Roi, etc.

### How to run
How to Run
Make sure you have:

`Fiji/ImageJ` with `TrackMate` installed.

Your `.tif `video files and optional `.txt `ROI files ready.

A config file as shown as the example in this file.

This script requires the TrackMate Java plugin, which is only available within the Fiji (ImageJ) environment. You can run this pipeline in two ways: 

##### Option 1: Run Inside Fiji (Recommended)
This is the most compatible method, especially if you're already using Fiji for spot detection and ROI generation.

Step-by-step:
1. Open Fiji.
2. Go to Plugins → Scripting → Script Editor.
3. Choose language: Python (Jython) from the top menu.
4. Copy-paste the full script (or use File → Open to load your .py script).
5. Click Run.

The output will be saved in the folder:
AnalysisSMolFiesta/ inside your TIFF directory.

Make sure:
- The script and config file are updated with the correct file paths.
- Your Fiji installation has TrackMate installed and updated (available via the Fiji update manager).

##### Option 2: Run via Python with Fiji’s Jython Environment
If you prefer to run the script from the command line:

Requirements:
- Python environment set up to use Fiji's ImageJ classes.
- Fiji's JAR files included in your CLASSPATH.
Use PyImageJ if running in a standalone Python environment.
But this option may need extra setup for Java integration and is generally more complex than using Fiji directly.

Then run the script:
```bash
python Trackmate.py
````

#### Notes
- Keep file naming consistent between .tif and .txt ROI files.
- Always verify segmentation match; the script will stop if ROI/TIFF counts mismatch.
- Output CSVs are sorted and reusable for downstream analyses (SMOL)


####  Configuration File Reference
This .ini-style config file defines how tracking is performed and where inputs/outputs are stored. Each section below controls a specific aspect of the pipeline:

| Section       | Key                                | Description                                                                                   | Example Value                                        |
| ------------- | ---------------------------------- | --------------------------------------------------------------------------------------------- | ---------------------------------------------------- |
| `[Paths]`     | `dir_tiff`                         | Folder containing input `.tif` videos to be analyzed                                          | `C:\...\NoMetadata`                                  |
|               | `dir_outline`                      | Folder with ROI outline `.txt` files, one per TIFF (if segmentation is enabled)               | `C:\...\segmented`                                   |
|               | `segmentation_enabled`             | Set to `True` to use ROI-based segmentation via outlines; otherwise, entire image is analyzed | `True` / `False`                                     |
| `[Tracking]`  | `linking_max_distance`             | Max distance (in pixels) to link spots across frames                                          | `20.0`                                               |
|               | `gap_closing_max_distance`         | Max distance allowed when closing temporal gaps                                               | `2.0`                                                |
|               | `splitting_max_distance`           | Max distance for splitting a spot into two tracks                                             | `2.0`                                                |
|               | `merging_max_distance`             | Max distance for merging two tracks into one                                                  | `2.0`                                                |
|               | `alternative_linking_cost_factor`  | Factor for alternative linking cost function in TrackMate                                     | `10.05`                                              |
| `[Detection]` | `threshold`                        | Intensity threshold for spot detection                                                        | `2.0`                                                |
|               | `detector_type`                    | Spot detector type (`Dog`, `LoG`, `Thresh`, etc.)                                             | `Dog`                                                |
| `[General]`   | `analysis_suffix`                  | Optional suffix added to output file names                                                    | `_experimentA` or leave empty                        |
| `[ROI]`       | `expansion_size`                   | Size (in pixels) to expand ROI borders (to avoid clipping spots at edges)                     | `0.0`                                                |
|               | `roi_ids`                          | Which ROI indices to use: `'all'` or comma-separated list (`1,2,3`)                           | `all`                                                |
| `[Settings]`  | `radius`                           | Expected average radius of spots (in pixels) for detection                                    | `3.0`                                                |
|               | `max_frame_gap`                    | Max number of frames allowed for linking a broken track (gap closing)                         | `10`                                                 |
|               | `number_threads`                   | Number of threads for parallel processing in TrackMate                                        | `16`                                                 |
| `[Filters]`   | `spot_filter1`                     | Spot filter in format: `FEATURE, threshold, keep_above=True/False`                            | `QUALITY, 3.0, True`                                 |
|               | `track_filter1`, `track_filter2` | Track-level filters (e.g., by number of spots, duration, etc.)                                | `NUMBER_SPOTS, 4, True`, `TRACK_DURATION, 4, True` |

- All paths must be absolute. Use double backslashes on Windows or raw strings (e.g., r"C:\path").

- Multiple filters can be added (e.g., spot_filter2, track_filter3, etc.).

- The order of filters in [Filters] can affect which tracks/spots survive preprocessing.

### Troubleshooting

| Issue / Error Message                                    | Possible Cause                                                                          | Solution                                                                                                  |
| -------------------------------------------------------- | --------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| `Timeout: Function track_process exceeded 60 seconds.`   | Processing a large image or complex ROI took too long.                                  | Debug the `track_process` function interactively in Fiji. Try using smaller ROIs or reducing image size.  |
| `Input error in TrackMate configuration`                 | Invalid detector/tracker settings or missing parameters.                                | Double-check `Trackmate_configFile.txt` for typos, missing sections, or unsupported settings.             |
| `Processing error in TrackMate`                          | TrackMate failed mid-process, possibly due to incompatible detector or corrupted image. | Test your image manually in Fiji using TrackMate GUI. Try switching detectors or re-exporting the image.  |
| `File mismatch, aborted`                                 | Number of TIFF images doesn't match the number of outline `.txt` ROI files.             | Ensure every TIFF has a corresponding ROI outline file. If not using segmentation, disable it in config.  |
| Fiji crashes without error                               | Memory overload from high-res images or too many ROIs.                                  | Close unused images in Fiji, reduce ROI size, or increase Fiji's memory via `Edit > Options > Memory`.    |
| No tracks are saved                                      | TrackMate found no valid tracks (e.g., all tracks filtered out).                        | Check spot quality, speed, and track length filters. Try relaxing thresholds in the config file.          ||
| `IndexError` on accessing ROI                            | Trying to access ROI that doesn’t exist (e.g., index too high).                         | Check that `roi_ids` values do not exceed the number of polygons in the outline file.                     |
| `AttributeError: 'NoneType' object has no attribute ...` | Fiji plugins (like TrackMate) or spot features missing or not loaded properly.          | Make sure you have the **latest TrackMate version** installed and Fiji is updated.                        |
| `No display window appears (viewer not rendering)`       | HyperStackDisplayer might be suppressed or running headless.                            | Viewer is optional. Comment out `displayer.render()` if running headless or use for debugging only.       |
| `PermissionError` or `FileNotFoundError`                 | Output directory not writable, or file paths incorrect.                                 | Check directory permissions, and ensure all paths in the config file exist. Use absolute paths if needed. |
| Track IDs or intensities look wrong                      | Spot detection might be too sensitive or noisy.                                         | Adjust the `threshold`, `radius`, or try different `detector_type` in the config.                         |
| Output CSVs contain nothing or zeros                | Spot quality too low or tracker didn’t detect valid movement.                           |	Try visualizing the detection in Fiji. Lower the detection threshold, adjust linking distances, or inspect raw TIFF + mask visually.              |
| `ModuleNotFoundError` or `ImportError`      | Required Java classes (e.g., `ij.*`, `fiji.plugin.trackmate.*`) are not accessible from Python. | Ensure you're using Fiji's **Jython scripting environment** or running via Fiji’s **Script Editor (Plugins > Scripting > Script Editor)** with language set to *Python (Jython)*. Make sure your Fiji is updated and compatible with the modules that are used in the script.|
| TrackMate or other Fiji plugins are missing | Your Fiji installation does not have TrackMate or supporting plugins installed.                 | Use **Fiji Updater** (`Help > Update Fiji`) and install `TrackMate`, `Bio-Formats`, and all recommended dependencies.                                                             |
| Script is run in CPython (standard Python)  | CPython can’t access Fiji’s Java libraries like `ImagePlus`.                                    | This script **must be run inside Fiji**, not in a standalone Python terminal. Launch via Fiji's Script Editor.                                                                    |
| `configparser.NoSectionError` or `NoOptionError` | Missing section or key in `Trackmate_configFile.txt`.  | Use the example config structure and ensure every required section (e.g., `[Tracking]`, `[Detection]`, `[Paths]`) is present. |
| `ValueError` when parsing config                 | Wrong value types (e.g., string where float expected). | Use proper types (e.g., `linking_max_distance = 5.0`, not `five`).                                                            |
| Tracks are saved, but look incorrect | Incorrect ROI, faulty segmentation, or high noise.  | Verify `.txt` outline files visually. Also check for proper channel selection and filter settings.           |
| Output files are not saved or are empty | Spot/tracking data not returned or overwritten. | Check log output and verify that tracks were detected. Confirm save paths are correct and filenames don’t clash. |
| `TrackMate window freezes or crashes`   | Too many ROIs or very large TIFF stacks.        | Try processing smaller batches, fewer ROIs, or increase Fiji memory (via `Edit > Options > Memory & Threads`).   |

