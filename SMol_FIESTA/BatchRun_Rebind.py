import os
import time
import argparse
import tomllib
from SMol_FIESTA import visualizer
from SMol_FIESTA import cell_info
from SMol_FIESTA import track_sorting
from SMol_FIESTA import bound_classification
from SMol_FIESTA import gaps_and_fixes
from SMol_FIESTA import rebind_analysis
from SMol_FIESTA import rebind_MSD
from SMol_FIESTA import rebind_fixed_particle
from natsort import natsorted
import numpy as np
from skimage import io as imgio

# Run all scripts in module on a provided config
# TODO: Modify command-centric workflow with new scripts
def run_scripts():
    parser = argparse.ArgumentParser(
        prog='BatchRun_Rebind',
        description='Run all scripts in module based on the config file provided. Perform track classification and rebinding analysis',
        epilog='')
    parser.add_argument('-c', '--config', default=None, type=str,
                        help='Specify the location of the config file in toml format.')
    parser.add_argument('-d', '--check-files', action='store_true',
                        help='Check the file structure provided in config, no action will be performed.')
    args = parser.parse_args()
    config_path = args.config
    print('Running as __main__ \n\t-> config_path: ' + str(config_path))

    # open config file
    if not config_path:
        print('Unknown config file location. Please pass in a config file via -c option.')
        print('Attempting to find default script-config.toml.')
        config_path = os.path.join(os.getcwd(), 'script-config.toml')
    try:
        with open(config_path, 'rb') as config_file:
            configs = tomllib.load(config_file)
    except FileNotFoundError:
        print('Config file not found or invalid.')
        return

    # only checks the file structure
    if args.check_files:
        print('\n\"check-files\" option is selected. Running checks on the file structure specified in config.')
        csv_path = configs['path']['csv_path']
        mask_path = configs['path']['mask_path']
        if configs['toggle']['run_fixed_particle']:
            comdet_path = configs['path']['comdet_path']
            comdet_files = natsorted(track_sorting.get_file_names_with_ext(comdet_path, 'csv'))
        masks = natsorted(track_sorting.get_file_names_with_ext(mask_path, 'png'))
        csv_sorted = track_sorting.csv_name_sort_suffix(csv_path, 'spots')
        csv_keys = natsorted(list(csv_sorted.keys()))

        # check file list length match
        print('\nCount\nTracks (csv):', len(csv_sorted.keys()),
              '\nMasks (png):', len(masks))
        if configs['toggle']['run_fixed_particle']:
            print('ComDet (csv):', len(comdet_files))
        print('\n\nMask-Track')
        if not len(csv_sorted.keys()) == len(masks):
            raise ValueError('Different number of Masks and Videos, you may have a video without tracks')

        # check number of cells in each mask
        for i in range(len(masks)):
            mask = np.swapaxes(imgio.imread(masks[i]), 0, 1)
            n_cell = np.max(mask)
            if configs['toggle']['run_fixed_particle']:
                print('\tComDet:', comdet_files[i])
            print('\tMask:', masks[i], '; number of cells:', n_cell)
            spots_video = track_sorting.index_format(natsorted(csv_sorted[csv_keys[i]]), n_cell)
            for k in range(len(spots_video)):
                if not spots_video[k] is None:
                    print('\t\t- Cell ' + str(k+1) + ':', spots_video[k])

        print('\nCheck complete. No errors found.')
        return

    # run full scripts
    else:
        track_sorting.main(config_path)
        print("")
        cell_info.main(config_path)
        print("")
        bound_classification.main(config_path)
        print("")

        if configs['toggle']['use_gap_fixed']:
            gaps_and_fixes.main(config_path)
            print("")

        rebind_analysis.main(config_path)
        print("")

        if configs['toggle']['run_fixed_particle']:
            rebind_fixed_particle.main(config_path)
            print("")

        if configs['toggle']['run_visualizer']:
            visualizer.main(config_path)
            print("")

        if configs['toggle']['run_MSD_calculations']:
            rebind_MSD.main(config_path)
            print("")
    return


if __name__ == '__main__':
    start_time = time.time()
    run_scripts()
    print("\nALL SCRIPTS: --- %s seconds ---" % (time.time() - start_time))