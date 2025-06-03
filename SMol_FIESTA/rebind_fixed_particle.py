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
    Either of:
        {csv_path}/{output_folder_name}/gaps-and-fixes_decisions.csv
        {csv_path}/{output_folder_name}/bound_decisions.csv
    Determined by parameter: {use_gap_fixed}

Output:
    {csv_path}/{output_folder_name}/rebind-fixed-particle_rebinding-events.csv
    {csv_path}/{output_folder_name}/rebind-fixed-particle_all-events.csv
    {csv_path}/{output_folder_name}/rebind-fixed-particle_assigned-spots.pkl

Parameters:
    allowed_dist_to_spot: allowed distance to the closest spot coordinate to be counted as colocalized
    conditional:
        use_gap_fixed: Use tracks processed by gaps_and_fixes.py, if False use tracks only processed by bound_classification.py instead.
"""

import numpy as np
import pandas as pd
from skimage import io as imgio
import time
import os
from natsort import natsorted
import logging
import pickle
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
    allowed_dist_to_spot = configs['rebind-fixed-particle'].get('allowed_dist_to_spot', 'None')
    use_gap_fixed = configs['toggle']['use_gap_fixed']

    # set default parameter when not specified
    allowed_spot_overlap = 0.0
    if isinstance(allowed_dist_to_spot, str): allowed_dist_to_spot = 0.0

    # file lists
    comdet_files = natsorted(get_file_names_with_ext(comdet_path, '.csv'))
    mask_files = natsorted(get_file_names_with_ext(mask_path, '.png'))

    rebind_df = pd.read_csv(str(os.path.join(output_path, 'rebind-strict-event.csv')))
    rebind_df = rebind_df.loc[:, ~rebind_df.columns.str.contains('^Unnamed')]

    tracks = pd.read_csv(str(os.path.join(output_path, ('gaps-and-fixes_decisions.csv' if use_gap_fixed else 'bound_decisions.csv'))))
    tracks = tracks.loc[:, ~tracks.columns.str.contains('^Unnamed')]
    headers = tracks[['Video #', 'Cell', 'Track']].to_numpy()
    tracks = slice_tracks_by_video(tracks, headers)
    headers = np.unique(headers, axis=0)

    if not len(comdet_files) == len(mask_files):
        raise ValueError('Different number of Masks and ComDet output, you may have an image without spot.')

    # Recorded output
    out_rebinds = []
    out_track_events = []
    out_spots = {}

    for i in range(len(mask_files)):
        # read files
        print_log('Processing', os.path.split(mask_files[i])[1], '->', os.path.split(comdet_files[i])[1])
        mask = np.swapaxes(imgio.imread(mask_files[i]), 0, 1)
        n_cell = np.max(mask)
        comdet = pd.read_csv(comdet_files[i], index_col=0)
        rebind_events = rebind_df[rebind_df['Video #'] == (i+1)] # fetch rebinding events with video number
        tracks_video = tracks[i+1] # fetch tracks with video number

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
                spots_cell[spot_ncell].append((len(spots_cell[spot_ncell]) + 1, spot_center, spot_area, spot_bounds, spot_cords))
            else:
                spots_cell[spot_ncell] = [(1, spot_center, spot_area, spot_bounds, spot_cords)]
        print_log('\t: Number of rogue spots (spots not in cell mask):', _)

        # filter spots and cells based on overlap in the same cell
        _ = [] # pending eliminated cells
        for cell in spots_cell:
            if len(spots_cell[cell]) < 2:
                continue
            spots = spots_cell[cell]
            kill = False
            for j in range(len(spots_cell[cell]) - 1):
                for k in range(len(spots_cell[cell])):
                    overlap = 0
                    for cords in spots[j][4]:
                        if np.any(np.all(cords == spots[k][4], axis=1)):
                            overlap += 1
                    if (float(overlap) / spots[j][2] > allowed_spot_overlap
                        or float(overlap) / spots[k][2] > allowed_spot_overlap):
                        _.append(cell)
                        kill = True
                        break
                if kill: break

        for cell in _:
            del spots_cell[cell]
        print_log('\t: Number of cells eliminated because of spots overlap:', len(_))

        # Assign rebinding events and tracks to cells
        if len(rebind_events) < 1:
            print_log('\t: SKIP: no strict rebinding event recorded for video #' + str(i+1))
            continue
        rebinds_cell = {}
        tracks_cell = {}
        for event in [rebind_events.iloc[i].to_dict() for i in range(len(rebind_events))]:
            cell = int(event['Cell'])
            if cell in rebinds_cell:
                rebinds_cell[cell].append(event)
            else:
                rebinds_cell[cell] = [event]
        for track_raw in tracks_video:
            track = track_raw[['Frame', 'x', 'y', 'Bound']].values
            cell = track_raw.iloc[0]['Cell']
            if cell in tracks_cell:
                tracks_cell[cell].append(track)
            else:
                tracks_cell[cell] = [track]

        # Process rebinding events, assign spot number to event endpoints
        _ = 0 # number of rebind endpoints assigned (each event has two endpoints)
        __ = 0 # number of rebind endpoints total
        for cell in rebinds_cell:
            for event in rebinds_cell[cell]:
                __ += 2
                s1 = pos_to_spot(event['x1'], event['y1'], spots_cell[cell], allowed_dist_to_spot, allowed_spot_overlap)
                s2 = pos_to_spot(event['x2'], event['y2'], spots_cell[cell], allowed_dist_to_spot, allowed_spot_overlap)
                if s1 != 0: _ += 1
                if s2 != 0: _ += 1
                event['From'] = s1
                event['To'] = s2
        print_log('\t: Assigned rebind endpoints:', _, 'Total', __)

        # Process track events, assign spot number to each event
        events_cell = {}
        _ = 0 # number of events assigned
        __ = 0 # number of events total
        for cell in tracks_cell:
            events = {}
            for track in tracks_cell[cell]:
                for event in split_track_to_behaviors(track):
                    event_avg = np.average(event, axis=0)
                    event_spot = pos_to_spot(event_avg[1], event_avg[2], spots_cell[cell], allowed_dist_to_spot, allowed_spot_overlap)
                    __ += 1
                    if event_spot != 0: _ += 1
                    if event_spot in events:
                        events[event_spot].append(event)
                    else:
                        events[event_spot] = [event]
            events_cell[cell] = events
        print_log('\t: Assigned track events:', _, 'Total', __)

        # Record rebind time and categorize each rebinding event
        rebind_times = {'same': [], 'diff': [], 'other': []}
        for cell in rebinds_cell:
            for event in rebinds_cell[cell]:
                if event['From'] == 0 or event['To'] == 0:
                    event['Spot'] = 'other'
                    rebind_times['other'].append(event['Time'])
                elif event['From'] == event['To']:
                    event['Spot'] = 'same'
                    rebind_times['same'].append(event['Time'])
                else:
                    event['Spot'] = 'diff'
                    rebind_times['diff'].append(event['Time'])
                out_rebinds.append(event)
        print_log('\tRebinding to the same spot:')
        if len(rebind_times['same']) != 0:
            print_log('\t-> ' + str(pd.Series(rebind_times['same']).describe()).replace('\n', '\t'))
        else:
            print_log('\t-> No event recorded.')
        print_log('\tRebinding to different spots:')
        if len(rebind_times['diff']) != 0:
            print_log('\t-> ' + str(pd.Series(rebind_times['diff']).describe()).replace('\n', '\t'))
        else:
            print_log('\t-> No event recorded.')
        print_log('\tRebinding (Others):')
        if len(rebind_times['other']) != 0:
            print_log('\t-> ' + str(pd.Series(rebind_times['other']).describe()).replace('\n', '\t'))
        else:
            print_log('\t-> No event recorded.')

        # Record bound time and categorize each track event
        bound_time = {'on': [], 'off': []} # on spot, off spot
        events_time = 0 # all event time
        for cell in events_cell:
            for spot in events_cell[cell]:
                for event in events_cell[cell][spot]:
                    behavior = int(event[0, 3])
                    event_avg = np.average(event, axis=0)
                    event_out = {'Video #': i+1, 'Cell': cell, 'Bound': behavior, 'Spot #': spot,
                                 'x_avg': event_avg[1], 'y_avg': event_avg[2], 'Time': event.shape[0]}
                    events_time += event.shape[0]
                    if behavior >= 2:
                        if spot == 0:
                            bound_time['off'].append(event.shape[0])
                        else:
                            bound_time['on'].append(event.shape[0])
                    out_track_events.append(event_out)
        print_log('\tBinding events colocalization with spot')
        print_log('\t-> on spot:', np.sum(bound_time['on']), 'frames')
        print_log('\t-> off spot:', np.sum(bound_time['off']), 'frames')
        print_log('\t-> total binding:', np.sum(bound_time['on']) + np.sum(bound_time['off']), 'frames')
        print_log('\t-> total events:', events_time, 'frames')
        print_log('\t-> proportion of on-spot binding / all binding:',
                  float(np.sum(bound_time['on'])) / (np.sum(bound_time['on']) + np.sum(bound_time['off'])))
        print_log('\t-> proportion of on-spot binding / all events:',
                  float(np.sum(bound_time['on'])) / events_time)

        # save spots in output
        out_spots[mask_files[i]] = spots_cell

    # summary stats
    print_log('\n[Analysis]')
    out_rebinds = pd.DataFrame(out_rebinds)
    out_track_events = pd.DataFrame(out_track_events)

    rebinds_same = out_rebinds.loc[out_rebinds['Spot'] == 'same']['Time'].to_numpy()
    rebinds_diff = out_rebinds.loc[out_rebinds['Spot'] == 'diff']['Time'].to_numpy()
    rebinds_other =  out_rebinds.loc[out_rebinds['Spot'] == 'other']['Time'].to_numpy()
    print_log('__________Rebind_________')
    print_log('Rebinding to the same spot:')
    if len(rebinds_same) != 0:
        print_log('-> ' + str(pd.Series(rebinds_same).describe()).replace('\n', '\t'))
    else:
        print_log('-> No event recorded.')
    print_log('Rebinding to different spots:')
    if len(rebinds_diff) != 0:
        print_log('-> ' + str(pd.Series(rebinds_diff).describe()).replace('\n', '\t'))
    else:
        print_log('-> No event recorded.')
    print_log('Rebinding (Others):')
    if len(rebinds_other) != 0:
        print_log('-> ' + str(pd.Series(rebinds_other).describe()).replace('\n', '\t'))
    else:
        print_log('-> No event recorded.')

    bound_events = out_track_events.loc[out_track_events['Bound'] >= 2]
    on_bound_time = np.sum(bound_events.loc[bound_events['Spot #'] > 0]['Time'].to_numpy())
    off_bound_time = np.sum(bound_events.loc[bound_events['Spot #'] == 0]['Time'].to_numpy())
    total_event_time = np.sum(out_track_events['Time'].to_numpy())
    print_log('\n____Binding_Colocalize___')
    print_log('on spot:', on_bound_time, 'frames')
    print_log('off spot:', off_bound_time, 'frames')
    print_log('total binding:', on_bound_time + off_bound_time, 'frames')
    print_log('total events:', total_event_time, 'frames')
    print_log('proportion of on-spot binding / all binding:', float(on_bound_time) / (on_bound_time + off_bound_time))
    print_log('proportion of on-spot binding / all events:', float(on_bound_time) / total_event_time)

    # output
    out_rebinds.to_csv(str(os.path.join(output_path, 'rebind-fixed-particle_rebinding-events.csv')))
    out_track_events.to_csv(str(os.path.join(output_path, 'rebind-fixed-particle_all-events.csv')))
    with open(str(os.path.join(output_path, 'rebind-fixed-particle_assigned-spots.pkl')), 'wb') as f:
        pickle.dump(out_spots, f, protocol=pickle.HIGHEST_PROTOCOL)
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
    for i in range(max_iterable):
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

def pos_to_spot(x:float, y:float, spots: list, dist_allowance:float = 0, allowed_overlap:float = 0):
    '''
    Assign the spot corresponding to the position given
    :param x: (float) x-coordinate
    :param y: (float) y-coordinate
    :param spots: (list) parsed spot information, (spot number, center, area, bounds, coordinates)
    :param dist_allowance: (float) distance closest to the spot coordinates that makes the position counted as in the spot
    :param allowed_overlap: (float) percentage area overlap allowed for particle spots in the same cell
    :return: (int) assigned spot number, 0 for background
    '''

    pix = np.array([x, y]).astype(int)

    # simplest case, don't need any calculations
    if allowed_overlap == 0 and dist_allowance == 0:
        for spot in spots:
            spot_coord = spot[4]
            if np.any(np.all(pix == spot_coord, axis=1)):
                return spot[0]
        return 0

    # go through each spot for matches
    pending = []
    for spot in spots:
        spot_coord = spot[4]
        distances = np.linalg.norm(spot_coord - np.array((x, y)), axis=1)
        min_dist = np.min(distances) # min dist to spot periphery
        if min_dist <= dist_allowance:
            pending.append(spot)
    if len(pending) == 0:
        return 0
    if len(pending) == 1:
        return pending[0][0]

    # if more than 1 spot is contested, evaluate based on distance to spot center
    min_dist = float('inf') # to center this time, might have position in the overlapping region
    winner = pending[0]
    for spot in pending:
        spot_dist = np.linalg.norm(np.array(spot[1])-np.array((x, y)))
        if spot_dist < min_dist:
            min_dist = spot_dist
            winner = spot
    return winner[0]

'''
================================================================================================================
TRACK PROCESSING
================================================================================================================
'''

# slice tracks based on headers
def slice_tracks_by_video(tracks, headers):
    indices = []
    save = np.array([-1, -1, -1])
    for i in range(headers.shape[0]):
        if not np.all(headers[i] == save):
            save = headers[i].copy()
            indices.append(i)
    indices.append(headers.shape[0])

    tracks_sliced = {}
    for i in range(len(indices) - 1):
        sliced = tracks.iloc[indices[i] : indices[i+1], :]
        video_n = sliced.iloc[0]['Video #']
        if video_n in tracks_sliced:
            tracks_sliced[video_n].append(sliced)
        else:
            tracks_sliced[video_n] = [sliced]
    return tracks_sliced


# split track by bound behavior, used for colocalization analysis
def split_track_to_behaviors(track:np.ndarray):
    indices = [] # indices to split the tracks, first frame of the next event
    behavior = track[0, 3] # first behavior
    for i in range(track.shape[0]):
        if track[i, 3] != behavior:
            behavior = track[i, 3]
            indices.append(i)
    return np.split(track, indices, axis=0)

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
    os.system('cls' if os.name == 'nt' else 'clear')


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