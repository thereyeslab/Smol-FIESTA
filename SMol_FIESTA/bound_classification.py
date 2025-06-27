"""
Runnable Script if run as __main__
Use the parsed tracks in track-sorting to determine track behavior by frame.
- Determine the relative behavior of spots given by their distances.
- Determine the behavior of the particle at each frame in the track based on spatial and temporal distances to
  neighboring frames.

Input:
    {csv_path}/{output_folder_name}/tracks.csv

Output:
    {csv_path}/{output_folder_name}/bound_decisions.csv

Parameters:
    distance_threshold: distance (pixels) for two spots to be considered within range for constrained diffusion.
    distance_threshold_strict: distance (pixels) for two spots to be considered within range for strict binding.
"""

import numpy as np
import pandas as pd
import time
import os
import logging
from tqdm import trange
import tomllib
import argparse

'''
================================================================================================================
CUSTOMIZABLE CRITERIA
================================================================================================================
'''
def determine_constrained(before:list, after:list):
    '''
        Determine given a list of binary criteria if the spot behavior is constrained diffusion.
        :param before: list(int[0,1]) criteria after the spots
        :param after: list(int[0,1]) criteria after the spots
        :return: int([0,1]) 1 if constrained diffusion, 0 if free diffusion
        '''
    if before[0] == 1 and after[0] ==1:
        return 1
    elif before[0] == 1 or after[0] == 1:
        return 1 if np.average(before[1:4]) > 0.75 or np.average(after[1:4]) > 0.75 else 0
    else: return 0

# Customizable: Criteria for strictly bound.
def determine_bound(before:list, after:list):
    '''
    Determine given a list of binary criteria if the spot behavior is bound.
    :param before: list(int[0,1]) criteria after the spots
    :param after: list(int[0,1]) criteria after the spots
    :return: int([0,1]) 1 if bound, 0 if diffusion
    '''
    decision = np.max([
        float(np.average(before[0:4])), float(np.average(after[0:4])),
        1.0 if np.average(np.concatenate([before[0:2], after[0:2]])) == 1.0 else 0.0,
        1.0 if np.average(np.concatenate([before[0:1], after[0:3]])) == 1.0 else 0.0,
        1.0 if np.average(np.concatenate([before[0:3], after[0:1]])) == 1.0 else 0.0,
        1.0 if np.average(before[1:4]) == 1.0 else 0.0,
        1.0 if np.average(after[1:4]) == 1.0 else 0.0,
    ])
    return 1 if decision >= 1 else 0

'''
================================================================================================================
MAIN
================================================================================================================
'''
def main(config_path:str = None):
    if not config_path:
        __location__ = os.path.realpath(
            os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, 'script-config.toml')
    with open(config_path, 'rb') as config_file:
        configs = tomllib.load(config_file)

    csv_path = configs['path']['csv_path']

    # Some parameters
    distance_threshold = configs['bound-classification']['distance_threshold']
    distance_threshold_strict = configs['bound-classification']['distance_threshold_strict']

    output_path = str(os.path.join(csv_path, configs['path']['output_folder_name']))
    if not os.path.isdir(output_path):
        raise ValueError('Directory do not exist, please run track_sorting.py first.')
    logging_setup(output_path, 'bound-classification')

    print_log('Reading from csv:', str(os.path.join(output_path, "Intermediates", 'tracks.csv')))
    tracks = pd.read_csv(str(os.path.join(output_path, "Intermediates", 'tracks.csv')))
    tracks = tracks.loc[:, ~tracks.columns.str.contains('^Unnamed')]

    headers = tracks[['Video #', 'Cell', 'Track']].to_numpy()
    tracks = slice_tracks(tracks, headers)

    track_results_threshold = []
    track_results_threshold_strict = []
    decisions_bound = []
    decisions_bound_strict = []
    decisions_bound_cdiffusion = []
    decisions_bound_overall = []
    for i in trange(len(tracks)):
        track = tracks[i]
        # Relaxed
        result_before, result_after = process_track(track, distance_threshold)
        track_results_threshold.append(result_before[result_before.columns[::-1]].join(result_after))
        decision_bound = []
        for j in range(len(result_before.index)):
            spot_before = result_before.iloc[j].to_numpy()
            spot_after = result_after.iloc[j].to_numpy()
            decision_bound.append(determine_constrained(spot_before, spot_after))
        decisions_bound.append(pd.Series(decision_bound))

        # Strict
        result_before, result_after = process_track(track, distance_threshold_strict)
        track_results_threshold_strict.append(result_before[result_before.columns[::-1]].join(result_after))
        decision_bound_strict = []
        record = [np.array([]), np.array([])]
        for j in range(len(result_before.index)):
            spot_before = result_before.iloc[j].to_numpy()
            spot_after = result_after.iloc[j].to_numpy()
            decision_bound_strict.append(determine_bound(spot_before, spot_after))
            record = [spot_before.copy(), spot_after.copy()]
        decisions_bound_strict.append(pd.Series(decision_bound_strict))

        # Constricted Diffusion
        decisions_bound_cdiffusion.append(pd.Series(np.logical_xor(decision_bound, decision_bound_strict).astype('int')))

        # Overall Bound 0, 1, 2
        decisions_bound_overall.append(pd.Series(np.add(decision_bound, decision_bound_strict)))

    output_relaxed = []
    output_strict = []
    output_both = []

    for i in trange(len(tracks)):
        track = tracks[i]
        result_relaxed = track_results_threshold[i]
        decision_relaxed = decisions_bound[i]
        output_relaxed.append(
            track[track.columns.tolist()[0:5] + track.columns.tolist()[8:]].reset_index(drop=True)
            .join(result_relaxed.assign(Bound=decision_relaxed).reset_index(drop=True))
        )
        result_strict = track_results_threshold_strict[i]
        decision_strict = decisions_bound_strict[i]
        output_strict.append(
            track[track.columns.tolist()[0:5] + track.columns.tolist()[8:]].reset_index(drop=True)
            .join(result_strict.assign(Bound=decision_strict).reset_index(drop=True))
        )
        decision_cdiffusion = decisions_bound_cdiffusion[i]
        decision_overall = decisions_bound_overall[i]
        output_both.append(
            track[track.columns.tolist()[:8]].reset_index(drop=True)
            .assign(RelaxedBound=decision_relaxed, StrictBound=decision_strict)
            .assign(ContrictedDiffusion=decision_cdiffusion)
            .assign(Bound=decision_overall)
        )
    pd.concat(output_both).to_csv(str(os.path.join(output_path, 'Intermediates', 'bound_decisions.csv')))
    return

'''
================================================================================================================
TRACKS
================================================================================================================
'''

# Transform distances values for each spot to binary decisions based on threshold
def process_track(track, threshold:float):
    dist_columns = [c for c in track.columns if 'Dist' in c]
    dist_range = int(len(dist_columns) / 2.0)
    dist_columns_before = dist_columns[:dist_range][::-1]
    dist_columns_after = dist_columns[dist_range:]

    dist_before = track[dist_columns_before].to_numpy()
    dist_after = track[dist_columns_after].to_numpy()

    result_before_columns = ['B' + str(i+1) for i in range(dist_range)]
    result_after_columns = ['A' + str(i+1) for i in range(dist_range)]

    result_before = pd.DataFrame((dist_before < threshold).astype('int'), columns=result_before_columns)
    result_after = pd.DataFrame((dist_after < threshold).astype('int'), columns=result_after_columns)
    return result_before, result_after

# slice the dataframe of tracks into individual track
def slice_tracks(tracks, headers):
    indices = []
    save = np.array([-1, -1, -1])
    for i in range(headers.shape[0]):
        if not np.all(headers[i] == save):
            save = headers[i].copy()
            indices.append(i)
    indices.append(headers.shape[0])

    tracks_sliced = []
    for i in range(len(indices) - 1):
        tracks_sliced.append(tracks.iloc[indices[i] : indices[i+1], :])
    return tracks_sliced


'''
================================================================================================================
START
================================================================================================================
'''


# Setup Logging
def logging_setup(path:str, script_name:str):
    log_file = str(os.path.join(path, 'logs', 'LOG_' + script_name + '.txt'))
    log_targets = [logging.FileHandler(log_file)]
    logging.basicConfig(format='%(message)s', level=logging.INFO, handlers=log_targets, force=True)
    logging.StreamHandler.terminator = ''
    open(log_file, 'w').close()
    os.system('cls' if os.name == 'nt' else 'clear')


# Modified print
def print_log(*args, end='\n'):
    print(' '.join([str(a) for a in args]), end=end)
    logging.info(' '.join([str(a) for a in args] + [end]))


# Start Script
if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(
        prog='bound-classification',
        description='Use the parsed tracks in track-sorting to determine track behavior by frame.',
        epilog='Prior scripts: track_sorting.py')
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__ \n\t-> config_path: ' + str(args.config))
    main(args.config)
    print_log("--- %s seconds ---" % (time.time() - start_time))