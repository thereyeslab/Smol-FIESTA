#!/bin/bash

# Activate conda environment (update the path if needed)
source ~/mambaforge/etc/profile.d/conda.sh
conda activate smol_env # put the name of your conda environment here (e.g., smol_env)

# Run the pipeline
python -m SMol_FIESTA.BatchRun_Rebind -c "/Users/masoomeshafiee/Desktop/test/pablo/SF2script-config.toml" #put the path to your config file here (e.g., "/Users/masoomeshafiee/Desktop/test/pablo/SF2script-config.toml")

# Keep terminal open after finishing
read -p "Press [Enter] to close this window..."
