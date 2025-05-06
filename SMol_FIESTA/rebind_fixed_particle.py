"""
Runnable Script if run as __main__
Incorporating ComDet fixed particle detection results to rebinding analysis
- Interpreting ComDet fixed particle results and assign particles to cell masks.
- Classify rebinding events based on their origin and destination particle.
- Isolate and calculate rebinding characteristics with new classification

Input:
    {comdet_path}/*_Cell_*_spotsAll.csv
    {csv_path}/{output_folder_name}/rebind-strict-event.csv

Output:
    TODO: annotate when done

Parameters:
    TODO: annotate when done
    allowed_spot_overlap: percentage area overlap allowed for particle spots in the same cell, skips cell if exceeds.
"""

import numpy as np
import pandas as pd
import time
import os
import logging
import csv
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

    # TODO: essentially two-part script (1. assign particles to each cell; 2. assign rebind events to particles)

    return

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