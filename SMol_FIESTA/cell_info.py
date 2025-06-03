"""
Runnable Script if run as __main__
Provide basic info for the cells in segmented masks.
- Determine numbers of cells, lengths of each cell and area of each cell for each mask.

Input:
    {mask_path}/*.png

Output:
    {csv_path}/{output_folder_name}/_cell-info.csv

Parameters:
    N/A
"""

import numpy as np
import pandas as pd
from skimage import io as imgio
import time
import os
from natsort import natsorted
import tomllib
import argparse

def main(config_path:str = None):
    if not config_path:
        __location__ = os.path.realpath(
            os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, 'script-config.toml')
    with open(config_path, 'rb') as config_file:
        configs = tomllib.load(config_file)

    csv_path = configs['path']['csv_path']
    mask_path = configs['path']['mask_path']
    output_folder_name = configs['path']['output_folder_name']

    masks = natsorted(get_file_names_with_ext(mask_path, 'png'))
    table = []
    columns = ['Mask #', 'Mask Name', '# Cells', 'Cell', 'Area', 'Length']
    for i in range(len(masks)):
        mask = np.swapaxes(imgio.imread(masks[i]), 0, 1)
        n_cell = np.max(mask)
        sizes = np.unique(mask, return_counts=True)
        cells = tabulate_cells(mask, n_cell)
        for j in range(n_cell):
            cell = np.array(cells[j])
            minx, maxx, miny, maxy = min(cell[:, 0]), max(cell[:, 0]), min(cell[:, 1]), max(cell[:, 1])
            l = np.sqrt(np.power(maxx-minx+1, 2) + np.power(maxy-miny+1, 2))
            table.append([i + 1, os.path.basename(masks[i]), n_cell, j + 1, sizes[1][j + 1], l])
    table = pd.DataFrame(table, columns=columns)
    table.to_csv(str(os.path.join(mask_path, '_cell-info.csv')))
    table.to_csv(str(os.path.join(csv_path, output_folder_name, '_cell-info.csv')))
    return

# Tabulating pixels for each cell, listing all pixel coordinates for each cell.
def tabulate_cells(mask, n_cell):
    result = [[] for i in range(n_cell)]
    for x in range(mask.shape[0]):
        for y in range(mask.shape[1]):
            try:
                cell = int(mask[x][y])
                result[cell - 1].append((x, y))
            except:
                continue
    return result

def get_file_names_with_ext(path: str, ext: str):
    flist = []
    for root, dirs, files in os.walk(path):
        for file in files:
            fname = file.split('.')
            if (fname[-1] == ext):
                flist.append(os.path.join(root, file))
    return flist

if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(
        prog='cell-info',
        description='Provide basic info for the cells in segmented masks.',
        epilog='')
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__ \n\t-> config_path: ' + str(args.config))
    main(args.config)
    print("--- %s seconds ---" % (time.time() - start_time))