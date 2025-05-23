# Troubleshooting
This guide provides solutions to common issues encountered while running the Rebinding Analysis Pipeline.

### Configuration & Setup

| Issue                                                           | Possible Cause                                                                             | Suggested Fix                                                                                                                                          |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **1. `FileNotFoundError: script-config.toml`**                  | - No config file was provided and `script-config.toml` is missing in the script directory. | Ensure the config file is named `script-config.toml` and placed in the same directory as the script, or use `-c <path/to/config>` to specify a path to your config file. |
| **2. `tomllib.TOMLDecodeError`**                                | - The `.toml` config file is improperly formatted.                                         |  Check for missing quotes, equals signs, or indentation issues. You can validate your TOML with [toml-lint](https://www.toml-lint.com/).             |
| **3. `ValueError: Different number of Masks and Videos`**       | - The number of `.csv` tracking files does not match the number of mask `.png` files.      |  Ensure each video has a corresponding mask file. Check naming consistency and file count.                                                           |                                                       |
| **4. `ModuleNotFoundError: No module named 'SMol_FIESTA'`**     | - The SMol_FIESTA or submodules are not installed or not in the Python path.                   | Ensure `SMol_FIESTA` is in your working directory or `PYTHONPATH`.          |
| **5. `ImportError` due to relative imports**                    | - Running the script from a different directory or in an incompatible environment.         | Always run the script from the top-level directory of the project. You can also add `.` or the project root to your `PYTHONPATH`.                   |
| **6. `KeyError: 'path'` or `'toggle'` in config**               | - Missing keys in the TOML config file.                                                    |  Double-check the TOML config has `[path]` and `[toggle]` sections with the required keys.                                                           |                                                 |
| **7. `PermissionError` writing to output files**                | - Lack of write permissions in the output directory.                                       | Run with appropriate permissions, or change output paths in the config to a writable location.                                                      |
| **8. Script silently returns early**                           | - `--check-files` flag is used, but user expects processing.                               | Make sure you're not passing `-d` or `--check-files` unless you only want to validate inputs.                                                        |
| **9. `FileNotFoundError` despite providing a path in the config file**                           | incorrect path file                               | Make sure all paths in config file reflect the current folder structure. If you moved or renamed files or folders after your last run, update config file accordingly.                                                        |
| 10.`ModuleNotFoundError: No module named 'pandas (or other libraries)'` | Dependencies not installed. | Run `pip install -r requirements.txt` to install all required packages. |


# Track Sorting Script Troubleshooting Guide

This table outlines common errors you may encounter when running the TrackMate track sorting script, along with their causes and recommended solutions.

| Error Message / Symptom                                                                 | Possible Cause                                                                 | Suggested Fix                                                                                      |
|------------------------------------------------------------------------------------------|------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `ValueError: Different number of Masks and Videos...`                                   | A mask exists without a corresponding TrackMate `spotsAll.csv` file         | Ensure each video has both a mask and correctly named `Cell_spotsAll.csv` files.                    |
| `ValueError: Track csv file name not formatted correctly: missing "Cell" in name.`      | Input filename does not follow expected format                              | Rename file to include `Cell_<index>`, e.g. `Video_Cell_3_spotsAll.csv`.                             |
| `IndexError: list index out of range occurred...`<br>`TypeError: cannot unpack...`      | Malformed or empty track passed to splitting/filtering function             | Ensure tracks contain valid spot data and are correctly extracted from CSV.                          |
| `Error loading CSV file <file>: could not convert string to float: '...'`               | CSV contains non-numeric or corrupted data                                  | Open the CSV and check for invalid formatting.                          |
| `# Cells in Video: 0` or pixel/IndexError during mask handling                          | Mask image is blank, corrupted, or mismatched in size                       | Verify that masks are properly generated and match the image size of corresponding tracks.          |
| `KeyError: 'path'` or similar when parsing TOML                                          | The TOML config is missing required keys                                    | Ensure `script-config.toml` has `[path]`, `[track-sorting]`, and `[toggle]` sections with all fields |                      |
| Script skips cells silently or log says `ALL SPOTS ELIMINATED`                          | Filtering parameters too strict; tracks fail QC                             | Maybe loosen thresholds like `allowed_gap_max`, `concurrent_max`, or lower `track_length_min`.            |
| Output folder contains no data or logs appear incomplete                                | Script may have exited early due to one of the above issues                 | Check the `track-sorting.log` for warnings or errors and fix root cause.                             |
| `ValueError: Track csv is for a cell index exceeding number of cells in the mask.`      | Track file references a cell index that is larger than the number of labeled cells in the mask | Ensure the mask has all expected labeled cells and that the `Cell_<index>` in the filename does not exceed this number. You can use tools like **Cellpose** to generate proper labeled masks. |


---

## Tips
- Mask should not be binary (0 or 1). Each cell should has a unique intiger value and the background should be 0. Use a segmentation tool like **Cellpose** to generate labeled masks where each cell has a unique ID.
- **Always check `track-sorting.log`** in the output directory for detailed processing info.
- **Validate inputs** before running the script: masks, CSVs, and config parameters.
- **Consistent naming** is crucial—both masks and CSVs must use the same base name + correct suffix.

# Cell_info.py
| **Issue**                                                 | **Possible Cause**                                                                                             | **Suggested Fix**                                                                                                                  |
| --------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------                                          |
| **Script crashes with `FileNotFoundError`**               | `mask_path` or `output_folder_name` directory doesn't exist.                                                   | Ensure all necessary folders exist and correct in the config file before running the script.                                                                      |
| **Output contains zero cells**                            | Mask is empty or incorrectly segmented (all pixels are zero).                                                  | Verify your mask files contain valid, non-zero labels. Use a tool like **Cellpose** to generate proper masks with unique cell IDs. |
| **Some masks produce fewer cells than expected**          | Mask labels are missing or not contiguous.                                                                     | Visualize the masks to ensure each cell has a unique positive integer ID.                                                          |
| **Script runs but produces incorrect area/length values** | Mask orientation might be incorrect due to axis swapping (`np.swapaxes`).                                      | Validate the axis order; if results seem off, consider toggling the `swapaxes` line (line 41) or checking your mask format.                  |


# bound_classification.py
| **Issue**                                                                     | **Why It Happens**                                                                                            | **Recommended Fix**                                                                                                                                                                                        |
| ----------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **1. Script crashes with `FileNotFoundError` for `tracks.csv`**               | The required input file (`tracks.csv`) is missing. This usually means `track_sorting.py` wasn’t run first.    | Ensure you run `track_sorting.py` before this script. Confirm that `{csv_path}/{output_folder_name}/tracks.csv` exists.                                                                                    |
| **2. Script crashes with `ValueError: Directory does not exist`**             | The path `{csv_path}/{output_folder_name}` from the config file doesn’t exist.                                | Check that the `csv_path` and `output_folder_name` are correct in your `script-config.toml`. Make sure the directory is created by earlier scripts.                                                        |
| **3. Script crashes with `KeyError` when reading config values**              | The expected keys like `csv_path` or `distance_threshold` are missing in the TOML config file.                | Double-check that your `script-config.toml` contains all required sections:<br>`[path]` → `csv_path`, `output_folder_name`<br>`[bound-classification]` → `distance_threshold`, `distance_threshold_strict` |
| **4. Script silently fails to log output**                                    | Logging is set up, but uses `os.system('cls')` which only works on Windows. On macOS/Linux it fails silently. | Replace `os.system('cls')` with cross-platform clearing logic or simply remove it. Also ensure logging file is writable.                                                                                   |
| **5. `ValueError` in `min()` or `np.average()` in custom criteria functions** | Happens if `process_track` returns all-zero rows, e.g., due to thresholds being too strict.                   | Change the distance thresholds                                                                                                                                                                |

# gaps_and_fixes.py
| **Issue** | **Possible Cause** | **Suggested Fix** |
|-----------|--------------------|-------------------|
| `KeyError: 'column_name'` | A column name in the `config.toml` does not exist in the CSV file. | Check if column names in the config match those in the CSV (`track_id`, `spot_id`, etc.). |
| `FileNotFoundError` | Input CSV or config path is incorrect. | Double-check file paths and extensions. Use absolute paths if needed. |
| Tracks aren’t split as expected | `max_frame_gap` or `max_distance_gap` may not be well tuned. | Adjust these values in the config.  |
| Too many tracks filtered out | `min_track_length` is too high. | Reduce `min_track_length`. |
| Output CSV is empty | All tracks were filtered out. | Check log messages and relax filtering settings. |
| Distance columns are all 0 or NaN | Tracks are too short or columns are missing. | Ensure tracks have ≥2 points and columns like `x`, `y`, `frame` exist. |
| `ValueError` during processing | Missing or malformed data. | Check for NaNs or unexpected values in your input CSV. |
| `ValueError: Directory does not exist, please run track_sorting.py first.` | The output directory has not been created. | Run `track_sorting.py` before this script to generate the expected folder structure. |

# Rebind_analysis.py
| Problem Description                                                                 | Possible Cause                                                                 | Suggested Solution                                                                                             |
|--------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
|`NameError: name '__file__' is not defined`                                          | Script is being run in an environment where `__file__` is undefined (e.g., Jupyter) |  run the script from the command line (`python rebind_analysis.py`).              |
| `FileNotFoundError` for `config.toml`                                                | The script can’t find `config.toml` in the expected location                   | Ensure `config.toml` is in the same directory as the script, or modify the script to accept a config path.   |
| `KeyError` when reading from config                                                  | Missing or misspelled keys in the `config.toml` file                           | Double-check that all required keys are present and spelled correctly in the TOML file.                      |
| Script finishes but no output is generated                                           | Output folder not created, or filters removed all data                         | Verify config filters. Check logs to confirm that tracks are read, filtered, and written properly.      |
| Distance computation looks incorrect or inconsistent                                 | Input CSV may have missing or malformed coordinates                            | Validate input data.                    |
| `ValueError: Directory does not exist, please run track_sorting.py first.` | The output directory has not been created. | Run `track_sorting.py` before this script to generate the expected folder structure. |

# Visualizer.py
| Problem Description                                                           | Possible Cause                                                                 | Suggested Solution                                                                                             |
|--------------------------------------------------------------------------------|--------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
 | `FileNotFoundError: 'script-config.toml'` or user-provided config not found   | Config file is missing or path is incorrect                                   | Ensure `script-config.toml` exists in the working directory or use `-c <path>` to specify the config path.   |
| `KeyError` or missing keys when parsing config                                | Required keys (e.g., `path`, `toggle`, `track-sorting`) are missing           | Review `script-config.toml` format and ensure all necessary sections and keys are included.                  |
 | No video output saved or `.tif` file not created                              | Output path does not exist, or no valid input data found                      | Ensure correct input files are placed in the `{csv_path}/{output_folder_name}` folder.                       |
 | Script says “Saved to…” but file is not in the folder                         | Output folder (`visualizer`) was not created or was mislocated                | Check if `visualizer/` was created inside the correct `{csv_path}/{output_folder_name}` directory.           |
 | `pandas.errors.ParserError` when loading CSV                                  | Malformed or corrupt input decision CSV file                                  | Verify CSV is complete and properly formatted; regenerate using `bound_classification.py` or `gaps_and_fixes.py`. |
 | `IndexError` or `ValueError` during `slice_tracks` or slicing logic           | CSV file is missing `Video #`, `Cell`, or `Track` columns                     | Ensure input decision CSVs have the expected columns and consistent data.                                    |
 | `FileNotFoundError` for mask or outline files                                 | Missing `.png` mask                | Make sure mask are correctly named and placed in the `mask_path` directory.                |
  | Colored tracks appear misaligned in the video                                 | Input mask and spot coordinates are from different resolutions/datasets       | Double-check that input masks match the scale and alignment of the tracking dataset.                         |
  | All spots appear the same color (or black)                                     | CSV may contain incorrect or undefined `Bound` column values                  | Inspect the CSV and verify that `Bound` column has valid values like 0, 1, 2.                                |
 | `AttributeError: 'NoneType' object has no attribute 'get'`                    | Missing `[visualizer]` section or incorrect TOML nesting                      | Confirm that `spot_diameter` and `tail_effect_enabled` are correctly nested under `[visualizer]` in the TOML.|
  | Script runs but generates empty black frames                                   | No valid tracks for the selected video/mask                                   | Make sure `Video #` in CSV matches index of the mask files (starting at 1).                                  |
  | Script crashes on certain machines due to dependencies                        | Missing packages like `skimage`, `tifffile`, `pandas`, or `natsort`           | Run `pip install -r requirements.txt` or manually install needed libraries.                                         |

# MSD_analysis.py
| **Issue**                                                                                       | **Possible Cause**                                                                | **Suggested Fix**                                                                                                  |
| ----------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| **1. `FileNotFoundError` when reading `bound_decisions.csv` or `gaps-and-fixes_decisions.csv`** | Output directory path or filename mismatch.                                       | Ensure `csv_path` in `script-config.toml` is correct and the expected file exists in the output folder.            |
| **2. `KeyError: 'Bound'` during processing**                                                    | The input CSV is missing the `Bound` column.                                      | Confirm that the input file (e.g., `bound_decisions.csv`) was generated correctly and contains the `Bound` column. |
| **3. No diffusion coefficients appear in output CSV**                                           | Subdivisions are skipped due to insufficient length (`< frames_per_subdivision`). | Check if your tracks have enough points. Consider lowering `frames_per_subdivision`.                               |
| **4. Output CSV contains `NaN` for Diffusion Coefficient or Alpha**                             | Subdivision has < 4 points or fitting failed.                                     | Ensure tracks are long enough, and try adjusting parameters.                             | Recheck cutoff/min thresholds (e.g., `cutoff_bound`, `min_fdiffusion`) in the script or TOML.                      |
| **5. `curve_fit` errors or bad Gaussian fits**                                                  | Poor histogram data or badly scaled bins.                                         | Make sure `Valid_Count` filtering isn't removing too many events. Try increasing `bin_scaling_factor`.             |
| **6. PDF is created but plots are blank**                                                       | All tracks were filtered out.                                                     | Confirm your filtering logic and thresholds aren't too strict. Use `print()` to debug event counts.                |
| **7. Weird log-scaled x-axis in plots**                                                         | Data may include `0` or near-zero values.                                         | Ensure all D values > 0. Check the logic that filters by `min_fdiffusion`, `min_bound`, etc.                       |
| **8. Warnings about invalid values in weights**                                                | Histogram is trying to plot zero-length slices.                                   | Check that each bound category has at least 1 valid event after filtering.                                         |



# Debugging Tips
- Use logging messages to trace the pipeline step-by-step.
- Add print() or pdb.set_trace() for in-depth inspection.
- Try running each module independently.

#  Still stuck?
Feel free to open an issue on GitHub with:
- The error message
- A description of what you were trying to do
- Your config file (if not sensitive)
- Screenshot or logs