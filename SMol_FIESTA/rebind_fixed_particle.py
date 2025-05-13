"""
Runnable Script if run as __main__
Incorporating ComDet fixed particle detection results to rebinding analysis
- Interpreting ComDet fixed particle results and assign particles to cell masks.
- Classify rebinding events based on their origin and destination particle.
- Isolate and calculate rebinding characteristics with new classification

Input:
    {comdet_path}/*_Results.csv (from comdet)
    {csv_path}/{output_folder_name}/rebind-strict-event.csv
    {mask_path}/*.png

Output:
    TODO: annotate when done

Parameters:
    TODO: annotate when done
    allowed_spot_overlap: percentage area overlap allowed for particle spots in the same cell, skips cell if exceeds.
"""

import numpy as np
import pandas as pd
from skimage import io as imgio
import time
import os
from natsort import natsorted
import logging
import shutil
import tomllib
import argparse

def main(config_path: str = None):
    if not config_path:
        __location__ = os.path.realpath(
            os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, 'script-config.toml')
    with open(config_path, 'rb') as config_file:
        configs = tomllib.load(config_file)

    # paths
    csv_path = configs['path']['csv_path']
    mask_path = configs['path']['mask_path']
    comdet_path = configs['path']['comdet_path']

    # logging setup
    output_path = str(os.path.join(csv_path, configs['path']['output_folder_name']))
    logging_setup(output_path, 'rebind-fixed')
    if not os.path.isdir(output_path):
        raise ValueError('Directory do not exist, please run track_sorting.py first.')

    # additional configs and parameters
    allowed_spot_overlap = configs['rebind-fixed-particle'].get('allowed_spot_overlap', 'None')

    # set default parameter when not specified
    if isinstance(allowed_spot_overlap, str): allowed_spot_overlap = 0.0

    # file lists
    comdet_files = natsorted(get_file_names_with_ext(comdet_path, '.csv'))
    mask_files = natsorted(get_file_names_with_ext(mask_path, '.png'))
    rebind_df = pd.read_csv(str(os.path.join(output_path, 'rebind-strict-event.csv')))
    rebind_df = rebind_df.loc[:, ~rebind_df.columns.str.contains('^Unnamed')]

    if not len(comdet_files) == len(mask_files):
        raise ValueError('Different number of Masks and ComDet output, you may have an image without spot.')

    for i in range(len(mask_files)):
        # read files
        print_log('Processing', os.path.split(mask_files[i])[1], '->', os.path.split(comdet_files[i])[1])
        mask = np.swapaxes(imgio.imread(mask_files[i]), 0, 1)
        n_cell = np.max(mask)
        comdet = pd.read_csv(comdet_files[i], index_col=0)
        rebind_events = rebind_df[rebind_df['Video #'] == (i+1)]

        # parse and assign comdet spots to cell
        fixed_spots = [comdet.iloc[i].to_dict() for i in range(len(comdet))]
        spots_cell = {}
        _ = 0 # counter for eliminated spots
        for spot in fixed_spots:
            spot_center = (spot['X_(px)'], spot['Y_(px)'])
            spot_ncell = mask[int(spot_center[0]), int(spot_center[1])]
            if spot_ncell == 0:
                _ += 1
                continue
            spot_bounds = (spot['xMin'], spot['yMin'], spot['xMax'], spot['yMax'])
            spot_area = spot['NArea']
            spot_cords = simulate_spot_coordinates(spot_center, spot_area, spot_bounds)
            if spot_ncell in spots_cell:
                spots_cell[spot_ncell].append((spot_center, spot_area, spot_bounds, spot_cords))
            else:
                spots_cell[spot_ncell] = [(spot_center, spot_area, spot_bounds, spot_cords)]
        print_log('\t: Number of rogue spots (spots not in cell mask):', _)

        # filter spots and cells based on overlap in the same cell
        _ = [] # pending eliminated cells
        for cell in spots_cell:
            if len(spots_cell[cell]) < 2:
                continue
            spots = spots_cell[cell]
            for j in range(len(spots_cell[cell]) - 1):
                for k in range(len(spots_cell[cell])):
                    overlap = 0
                    for cords in spots[j][3]:
                        if np.any(np.all(cords == spots[k][3], axis=1)):
                            overlap += 1
                    if (float(overlap) / spots[j][1] > allowed_spot_overlap
                        or float(overlap) / spots[k][1] > allowed_spot_overlap):
                        _.append(cell)
        for cell in _:
            del spots_cell[cell]
        print_log('\t: Number of cells eliminated because of spots overlap:', len(_))



    return


'''
================================================================================================================
FIXED-SPOTS PROCESSING
================================================================================================================
'''

def simulate_spot_coordinates(center: tuple, area: int, bounds: tuple, max_iterable:int = 128):
    '''
    Simulate coordinates for each comdet-detected spots, +x, +y, -x, -y direction order
    :param center: (float) center of detection -> only use int parts (x, y)
    :param area: (int) area of the spot
    :param bounds: (tuple(int)) xMin, yMin, xMax, yMax
    :param max_iterable: maximum number of iterable coordinates the code will use for generating direction steps
    :return: (list(int)) list of coordinates for the simulated spot
    '''

    n_pix = area - 1
    directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    step = 0
    pix = np.array(center).astype(int)
    coords = [pix.copy()]
    min_bounds = [bounds[0], bounds[1]]
    max_bounds = [bounds[2], bounds[3]]

    # generate direction indices for a spiral out expansion of coordinates
    dirindex = []
    for i in range(256):
        if i < 2:
            dirindex.append(i)
            continue
        i1 = int(i / 2.0) + 1
        dirindex += [i % len(directions)] * i1

    while n_pix > 0:
        pix = pix + directions[dirindex[step]]
        step += 1

        if ((pix >= min_bounds) & (pix <= max_bounds)).all():
            coords.append(pix.copy())
            n_pix -= 1
    return coords
'''
================================================================================================================
FILE HANDLING
================================================================================================================
'''

# get all files under path with extension
def get_file_names_with_ext(path: str, ext: str):
    flist = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if os.path.splitext(file)[1] == ext:
                flist.append(str(os.path.join(root, file)))
    return flist

'''
================================================================================================================
START
================================================================================================================
'''


# Setup Logging
def logging_setup(path:str, script_name:str):
    log_file = str(os.path.join(path, 'LOG_' + script_name + '.txt'))
    log_targets = [logging.FileHandler(log_file)]
    logging.basicConfig(format='%(message)s', level=logging.INFO, handlers=log_targets)
    logging.StreamHandler.terminator = ''
    open(log_file, 'w').close()
    os.system('cls')


# Modified print
def print_log(*args, end='\n'):
    print(' '.join([str(a) for a in args]), end=end)
    logging.info(' '.join([str(a) for a in args] + [end]))


# Start Script
if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(
        prog='rebind-analysis',
        description='Performs analysis on bound-classified tracks.',
        epilog='Prior scripts: bound_classification.py, gaps_and_fixes.py')
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__ \n\t-> config_path: ' + str(args.config))
    main(args.config)
    print_log("--- %s seconds ---" % (time.time() - start_time))