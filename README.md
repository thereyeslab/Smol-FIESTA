# SMol_FIESTA

A modular and reproducible microscopy video analysis pipeline to extract and quantify the binding behavoir of proteins in the single-molecule flourescnet microscopy data using Python. 
#### Main logic

This pipeline processes single-particle tracking data obtained from microscopy videos (e.g., using [TrackMate](https://www.sciencedirect.com/science/article/pii/S1046202316303346))



1. Tracking Initialization (Fiji TrackMate Script)
Use the included `TrackMateScriptv2.py` to automate batch tracking in Fiji via command line. You can extract and export tracking data using segmentation masks and configuration from a .toml file.

2. Threshold Calculation for Step Size Filtering
Run `Thresh_calculator.py` to compute a noise-adjusted diffusion threshold based on pixel size, diffusion coefficient, and vibration noise. This threshold is used later to define diffusive vs. bound behavior.

3. Track Processing and Binding Classification
The main analysis pipeline (BatchRun_Rebind.py) processes all tracks to:

- Filter short/noisy tracks

- For each spot in a molecule’s trajectory, it calculates the distance to neighboring spots in both past and future frames. Based on configurable thresholds and two classification criteria (strict and relaxed), it assigns each spot in a trajectory a behavioral state:

    - Bound

    - Constrained Diffusion

    - Diffusive

- Fill gaps and interpolate short missing segments: After the initial classification, additional filtering steps refine each track’s overall behavior by analyzing the type and duration of surrounding events.

- Quantify rebinding events: The pipeline then detects possible rebinding events within each trajectory and quantifies:
    - Average binding and rebinding durations
    - Frequency of state transitions
     Proportions of binding states and rebinding events at the population level

4. Kinetic Modeling via Two-Exponential Fit
After classification, `TwoExponential.py` fits the dwell-time distribution to either a single or two-exponential model. This provides kinetic parameters for each molecular state, including:

- Mean lifetimes (τ) (aka bound time or search time, depending on the inputs)

- Decay rates (λ)

- Transition probabilities

- State assignment confidence per event

For plots, module details, and examples, see the `docs/` and `example_data/` directories.


**Inputs**
1. Config File (.toml)
Specifies paths, parameters, and module toggles ( check the `config.md` for details)

2. TrackMate spots files
One CSV per video/cell. Required naming: *_Cell_*_spotsAll.csv
Typical columns include:
    | Column Name   | Description                              |
    |---------------|------------------------------------------|
    | TRACK_ID      | ID of the track     |
    | POSITION_T    | Time point (frame number).               |
    | POSITION_X    | X coordinate (in pixels).    |
    | POSITION_Y    | Y coordinate (in pixels).    |
    | INTENSITY     | Spot intensity             |

3. Segmentation Masks (.png)
Grayscale images with unique integer labels per cell (0 = background).
Recommended tool: [Cellpose](https://github.com/MouseLand/cellpose)

**Outputs**
See `docs/input_output.md` and the `example_data/` folder for full output samples. Key outputs:
- Threshold printout: From Thresh_calculator.py for filtering setup (Ancillary script)
- bound_decisions.csv — Spot-level binding state for each frame
- rebind-Events.csv - Duration, time and types of all events (binding, rebinding, constrained diffusion, fast diffusion, search time ) within each track
- RESULT_rebind.csv — Summary statistics of rebinding
- RESULT_rebind.txt - Summary statistics of binding behavior in more details
- 2Exponential_test.pdf: Fitted exponential models of the event durations and histograms
- Exp1_probabilities.csv: Posterior state probabilities per event
- transition_probabilities.csv: Transition matrix between states
- Diffusion_Coefficient_Calculation.csv — Track-level diffusion coefficients
- Diffusion_Coefficient_Plots.pdf — Plots of diffusion coefficient distributions
- Reconstructed videos — Simulated videos with color-coded spot behaviors

## Features

- Configurable and modular
- Track-to-cell assignment from segmentation masks
- Automated TrackMate ROI-based tracking via Fiji
- Step-size threshold computation with vibration noise support, to be used for classifying the molecules binding states. 
- Bound state classification based on two different cirteira (relaxed and strict)
- Gap filing and interpolation.
- Extra filtering steps for baceteria images. 
- Rebinding detection and analysis of transition between states (outputting transition matrix)
- Calculation of molecule's Search time for binding
- Molecule Rebinding and binding event analysis
- MSD-based diffusion coefficient estimation (Mean Squared Displacement)
- Two-exponential model fitting with posterior assignment to identify the mixed behaviors
- Tabular + visual outputs
- Statistical summaries
- Real-time logging to easily track errors and monitor pipeline progress
- Metadata saving for reproducibility

## My Contributions

As part of the **Reyes Lab** development team, I worked on both **code development** and **project infrastructure**, with the goal of making the pipeline more robust, reproducible, and user-friendly. My main contributions include:

- **Code development & refactoring** → contributed to the development of new functionality, reorganized existing modules, improved readability, and applied clean code principles. Suggested and implemented improvements to the overall workflow logic (e.g., feature extraction, motion classification, pipeline orchestration).  
- **Deep code review ** → went through all scripts in detail to fully map the workflow, verify correctness, and document the design logic at the module level.  
- **Documentation**  
  - Authored detailed module-level documentation explaining how each component of the pipeline works and how to use it.  
  - Created a dedicated document for **input/output data flow**, describing file formats, column definitions, and expected results.  
  - Wrote a **troubleshooting guide**, anticipating common user issues (especially for non-technical users) and presenting solutions in a structured table format per module.  
  - Authored a top-level **README.md** with setup instructions, usage examples, and overall guidance for new users.  
- **Cross-platform usability** → tested the pipeline on macOS and Windows, and created launcher scripts.
- **Testing & debugging** → ran the pipeline both end-to-end and module-by-module with example datasets to validate results.
- **User experience improvements** → improved how outputs and intermediate results are represented, making the analysis easier to interpret and more consistent across runs.  
- **Future development** → compiled a list of potential issues, enhancement ideas, and improvement opportunities to inform future versions of the project.  

Overall, my work combined **code contributions, refactoring, and documentation** to strengthen the pipeline’s usability and maintainability for both technical and non-technical users.


## Installation

**1. Download or clone the repository:**
``` bash
git clone https://github.com/thereyeslab/Smol-FIESTA.git
cd SMol-FIESTA
```
**2. Create and activate a conda environment:**
``` bash
conda create -n smol_env python=3.11
conda activate smol_env
```
**3. Install dependencies:**
``` bash
pip install -r requirements.txt
```
**4. Install as a package (Optional)**
```bash
pip install -e .
```
## Usage
1. Modify `script-config.toml` to match your data paths and settings. Consult the `config.md` for detail explanation. 
2. Activate your environment if not already. 
``` bash
conda activate smol_env
```
3. Run threshold calculator (optional but recommended). You can use this to get a threshold for the step size used for binding behavior classification. See thresh_calculator.md for how to use this script. 
4. Run TrackMateScriptv2 in batch mode in Fiji to get your tracks and spots files compatible to the pipeline. See TrackMateScriptV2.md for setup. You can skip this if you have spots.csv and tracks.csv files compatible with the pipeline. See input_outputs.md.
5. Run main pipeline: There are three different ways to run the pipeline: 

- A. Run via command line:
``` bash
python -m SMol_FIESTA.BatchRun_Rebind -c /path/to/your/script-config.toml
```
- B. Run the package ( If already installed with `pip install -e .` ):
``` bash
SMF -c /path/to/your/script-config.toml
```
- C. By double-clicking a launcher script
    - On Windows: Use run_smol.bat
    - On macOS: Use run_smol.command

After editing the script to insert your config file path and environment name:
```bash
chmod +x run_smol.command  # ONLY once on macOS. `.bat` file in Windows is already double-clickable
```
You can also create a desktop shortcut for ease of use.

##### Advanced Usage

Run a specific module:
``` bash
python -m SMol_FIESTA.Module_Name -c /path/to/your/script-config.toml
```
Check files only (no processing):
``` bash
python -m SMol_FIESTA.BatchRun_Rebind -c /path/to/your/script-config.toml --check-files
```

## Project Structure
``` r
SMol_FIESTA/
├── SMol_FIESTA/              # Main source code
│   ├── BatchRun_Rebind.py    # Main runner
│   ├── track_sorting.py
│   ├── cell_info.py
│   ├── bound_classification.py
│   ├── gaps_and_fixes.py
|   ├── rebind_analysis.py
|   ├── rebind_MSD.py
|   ├── TwoExponential.py
├── visualizer.py
| 
├── Ancillary scripts 
|   ├── Thresh_calculator.py
|   ├── TrackmateScriptv2.py
|   ├── thresh_calculator.md
|   ├── TrackMateScriptV2.md
|   ├── Trackmate_configFile.txt
|
├── docs                    # Helpfull documentation
├── requirements.txt
├── setup.py
├── run_smol.bat           # For Windows
├── run_smol.command       # For Mac OS
├── README.md
├── pyproject.toml
└── script-config.toml        # Example config file

```




## Example Usage

See the example_data/ directory for a working example with input files, config, and output results.
For common issues and solutions, check `docs/troubleshooting.md`.

## Configuration File
The .toml config defines file paths, thresholds, toggles, and options for each module.
It is copied to the output folder for reproducibility.
See `docs/config.md` for full documentation. 

## Acknowledgment

This work is based on the original SMol-FIESTA repository developed by the Reyes Lab.

## Contributors
- Jose Pablo Rascon Perez
- Steven Hu
- Masoumeh Shafiei

Contributions are welcome! Feel free to open issues or submit pull requests.

## Citation
If you use this tool for your research, please cite:
Reyes-Lamothe Lab, McGill University
[Author names, DOI, publication] # to be done

