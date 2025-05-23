# SMol_FIESTA

A modular and reproducible microscopy video analysis pipeline to extract and quantify the binding behavoir of proteins in the single-molecule flourescnet microscopy data using Python. 
#### Main logic

This pipeline processes single-particle tracking data obtained from microscopy videos (e.g., using [TrackMate](https://www.sciencedirect.com/science/article/pii/S1046202316303346))

For each spot in a molecule’s trajectory, it calculates the distance to neighboring spots in both past and future frames. Based on configurable thresholds and two classification criteria (strict and relaxed), each spot is assigned to one of the following states:
- Bound
- Constrained diffusion
- Diffusive

After the initial classification, additional filtering steps refine each track’s overall behavior by analyzing the type and duration of surrounding events. The pipeline then detects possible rebinding events within each trajectory and quantifies:
- Average binding and rebinding durations
- Frequency of state transitions
- Proportions of binding states and rebinding events at the population level

For diagrams, module details, and examples, see the `docs/` and `example_data/` directories.


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
# Decide later
See `docs/input_output.md` and the `example_data/` folder for full output samples. Key outputs:
- bound_decisions.csv — Spot-level binding state for each frame
- rebind-strict-boundtime.csv — Duration of binding events
- rebind-strict-rebindingtime.csv — Duration between binding events
- rebind-AllDiffusion-time.csv — Duration of all diffusion events
- RESULT_rebind.csv — Summary statistics of rebinding
- Diffusion_Coefficient_Calculation.csv — Track-level diffusion coefficients
- Diffusion_Coefficient_Plots.pdf — Plots of diffusion coefficient distributions
- Reconstructed videos — Simulated videos with color-coded spot behaviors

## Features

- Configurable and modular
- Track-to-cell assignment from segmentation masks
- Bound state classification based on two different cirteira (relaxed and strict)
- Gap filing and interpolation.
- Extra filtering steps for baceteria images. 
- Molecule Rebinding and binding event analysis
- Compute and plot MSD distribution (Mean Squared Displacement)
- Tabular + visual outputs
- Statistical summaries
- Real-time logging to easily track errors and monitor pipeline progress
- Metadata saving for reproducibility

## Installation

**1. Download or clone the repository:**
``` bash
git clone https://github.com/thereyeslab/Smol-FIESTA.git
cd SMol_FIESTA
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
3. There are three different ways to run the pipeline: 

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
|   ├── visualizer.py
|
├── docs                    # Helpfull documentation
├── requirements.txt
├── setup.py
├── srun_smol.bat           # For Windows
├── srun_smol.command       # For Mac OS
├── README.md
└── script-config.toml        # Example config file

```




## Example Usage

See the example_data/ directory for a working example with input files, config, and output results.
For common issues and solutions, check `docs/troubleshooting.md`.

## Configuration File
The .toml config defines file paths, thresholds, toggles, and options for each module.
It is copied to the output folder for reproducibility.
See `docs/config.md` for full documentation. 

## Contributors
- Jose Pablo Rascon Perez
- Steven Hu
- Masoumeh Shafiei

Contributions are welcome! Feel free to open issues or submit pull requests.

## Citation
If you use this tool for your research, please cite:
Reyes-Lamothe Lab, McGill University
[Author names, DOI, publication] # to be done

