import os
import time
import argparse
import tomllib
from SearchAnalysis import visualizer
from SearchAnalysis import cell_info
from SearchAnalysis import track_sorting
from SearchAnalysis import bound_classification
from SearchAnalysis import gaps_and_fixes
from SearchAnalysis import rebind_analysis
from SearchAnalysis import rebind_MSD
from natsort import natsorted
import numpy as np
from skimage import io as imgio

# Run all scripts in module on a provided config
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
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, 'script-config.toml')
    with open(config_path, 'rb') as config_file:
        configs = tomllib.load(config_file)

    # only checks the file structure
    if args.check_files:
        print('\n\"check-files\" option is selected. Running checks on the file structure specified in config.')
        csv_path = configs['path']['csv_path']
        mask_path = configs['path']['mask_path']
        masks = natsorted(track_sorting.get_file_names_with_ext(mask_path, 'png'))
        csv_sorted = track_sorting.csv_name_sort_suffix(csv_path, 'spotsAll')
        csv_keys = natsorted(list(csv_sorted.keys()))

        # check file list length match
        print('\nCount\nTracks (csv):', len(csv_sorted.keys()), '; Masks (png):', len(masks), '\n\nMask-Track')
        if not len(csv_sorted.keys()) == len(masks):
            raise ValueError('Different number of Masks and Videos, you may have a video without tracks')

        # check number of cells in each mask
        for i in range(len(masks)):
            mask = np.swapaxes(imgio.imread(masks[i]), 0, 1)
            n_cell = np.max(mask)
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