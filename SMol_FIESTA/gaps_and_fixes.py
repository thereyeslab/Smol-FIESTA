"""
Runnable Script if run as __main__
Further refining track behavior interpretation by interpolating the gap behavior.
- Determine gap behavior by temporal neighboring events, gap length and event lengths.
- Introduce conditional filtering by overall binding event probability.

Input:
    {csv_path}/{output_folder_name}/bound_decisions.csv

Output:
    {csv_path}/{output_folder_name}/gaps-and-fixes_decisions.csv


    IF filter_by_binding_prop:
        {csv_path}/{output_folder_name}/filtered_passed_spotsAll.csv
        {csv_path}/{output_folder_name}/filtered_failed_spotsAll.csv

Parameters:
    min_time_strict: min length (frames) of strict binding event.
    min_time_constrained: min length (frames) of constrained diffusion event.
    min_time_diffusion: min length (frames) of free diffusion event.
    max_bound_gapFill: max length (frames) of gap that can be interpreted as strict binding.
    min_prop_binding: filter condition after gap-fixes, eliminate free particle tracks.
    max_prop_binding: filter condition after gap-fixes, eliminate permanently-bound tracks.

    conditional:
        allowed_gapTotal: total gap allowed to be filled per track
        Prop_check_over
        Ratio_turnover_dtf
        Ratio_turnover_all
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

def main(config_path:str = None):
    if not config_path:
        __location__ = os.path.realpath(
            os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, 'script-config.toml')
    with open(config_path, 'rb') as config_file:
        configs = tomllib.load(config_file)

    csv_path = configs['path']['csv_path']

    # Some parameters
    filter_by_binding_prop = configs['toggle'].get('filter_by_binding_prop', 'None')
    min_time_strict = configs['gaps-and-fixes'].get('min_time_strict', 'None')
    min_time_constrained = configs['gaps-and-fixes'].get('min_time_constrained', 'None')
    min_time_diffusion = configs['gaps-and-fixes'].get('min_time_diffusion', 'None')
    max_bound_gapFill = configs['gaps-and-fixes'].get('max_bound_gapFill', 'None')
    min_prop_binding = configs['gaps-and-fixes'].get('min_prop_binding', 'None')
    max_prop_binding = configs['gaps-and-fixes'].get('max_prop_binding', 'None')
    allowed_gapTotal = configs['gaps-and-fixes'].get('allowed_gapTotal', 'None')


    # turning off unused filters

    if isinstance(filter_by_binding_prop, str): filter_by_binding_prop = False
    if isinstance(min_time_strict, str): min_time_strict = 0
    if isinstance(min_time_constrained, str): min_time_constrained = 0
    if isinstance(min_time_diffusion, str): min_time_diffusion = 0
    if isinstance(max_bound_gapFill, str): max_bound_gapFill = 0
    if isinstance(min_prop_binding, str): min_prop_binding = 0
    if isinstance(max_prop_binding, str): max_prop_binding = 2.0
    if isinstance(allowed_gapTotal, str): allowed_gapTotal = 1000000

    output_path = str(os.path.join(csv_path, configs['path']['output_folder_name']))
    if not os.path.isdir(output_path):
        raise ValueError('Directory do not exist, please run track_sorting.py first.')
    logging_setup(output_path, 'gaps-and-fixes')


    print_log('Reading from csv:', str(os.path.join(output_path, 'Intermediates', 'bound_decisions.csv')))
    tracks = pd.read_csv(str(os.path.join(output_path, 'intermediates', 'bound_decisions.csv')))
    tracks = tracks.loc[:, ~tracks.columns.str.contains('^Unnamed')]

    headers = tracks[['Video #', 'Cell', 'Track']].to_numpy()
    tracks = slice_tracks(tracks, headers)

    counts_gap = 0
    counts_event = 0
    counts_operations = 0
    counts_filtered = 0
    frames_filtered = 0
    frames_total = 0
    counts_filtered_onlyBound = 0
    output_tracks = []
    frames_filtered_onlyBound =0

    for i in range(len(tracks)):
        track = tracks[i]
        header = track[['Video #', 'Cell', 'Track']].to_numpy()[0]
        pos = track[['x', 'y']].to_numpy()
        vname = track.iloc[0]['Video Name']
        track = track.loc[:, ~track.columns.str.startswith(('Video Name', 'x', 'y'))]
        track = list(track.to_numpy())
        print_log('Fixing:', 'Video', header[0], 'Cell', header[1], 'Track', header[2])

        print_log('\t-> Track Duration:', len(pos), 'frames.')
        save_track = len(pos)
        # Fill Gaps
        _, track, pos = process_gaps(track, pos, lambda l, r, dur: 1 if l == 2 and r == 2 and dur > max_bound_gapFill else min(l, r))
        print_log('\t-> Gap:', _, 'filled')
        if _ > allowed_gapTotal:
            print_log('\t\t: * Gap counts within track exceeds threshold: Aborting Track')
            continue
        counts_gap += _
        # Separate Events
        _, events, dtftransition = event_separation(track)
        #if dtftransition == 0:
        #    continue
        print_log('\t-> Events:', _, 'found')
        found_events = _
        counts_event += _
        print_log('\t-> EventRatio All transitions:', found_events/save_track)
        print_log('\t-> EventRatio Diffusion and Bound :', dtftransition/save_track)

        __ = 0
        # Pass 3: CD -> FD, Merge
        _, events1 = pass_events(events, 1, lambda l,r: 2 if l == 2 and r == 2 else 1, min_time_constrained)
        print_log('\t-> Pass 3:', _, 'events relabeled.')
        __ += _

        # Pass 1: FD -> CD, Merge
        _, events2 = pass_events(events1, 0, lambda l,r: 2 if l == 2 and r == 2 else 1, min_time_diffusion)
        print_log('\t-> Pass 1:', _, 'events relabeled.')
        __ += _

        # Pass 2: SB -> CD, Merge
        _, events3 = pass_events(events2, 2, 1, min_time_strict)
        print_log('\t-> Pass 2:', _, 'events relabeled.')
        __ += _



        print_log('\t-> Pass Complete:', __, 'relabeling performed,', len(events3), 'events left.')
        counts_operations += __

        # Collect Events -> Reconstruct Track
        track1 = events_to_track(events1)
        track2 = events_to_track(events2)
        track3 = events_to_track(events3)

        # Filter by binding proportion, weighted
        if filter_by_binding_prop:
            prop_binding = prop_counting(track3, (1.0, 1.0))  # weights (constricted, strict), for min  #track1 w0
            prop_binding_only_strict = prop_counting(track3, (0, 1.0), False)  # weights for max
            print_log('\t-> Filter (MIN BINDING PROPORTION):', prop_binding, 'non-diffusing')
            print_log('\t-> Filter (MAX BINDING PROPORTION):', prop_binding_only_strict, 'strictly bound')

            frames_total += len(track1)

            if(prop_binding < min_prop_binding):
                print_log('\t\t: FAIL by minimum proportion binding')
                frames_filtered += len(track1)
                counts_filtered += 1

                continue
            elif(prop_binding_only_strict > max_prop_binding):
                print_log('\t\t: FAIL by maximum proportion binding')
                frames_filtered_onlyBound += len(track1)
                counts_filtered_onlyBound += 1
                continue
            else:
                print_log('\t\t: PASS')


        trackdf = pd.DataFrame(np.array(track3)[:, :5], columns=['Video #', 'Cell', 'Track', 'Frame', 'Intensity'])
        trackdf = (trackdf.assign(GapFixed=np.array(track)[:, 8],
                                  Pass1=np.array(track1)[:, 8], Pass2=np.array(track2)[:, 8], Pass3=np.array(track3)[:, 8],
                                 Bound=np.array(track3)[:, 8], isGap=np.array(track3)[:, 9],
                                 Name=np.array([vname]*len(track3)))
                   .join(pd.DataFrame(np.array(pos), columns=['x', 'y']))
                   .rename(columns={'Name':'Video Name'}))
        trackdf = trackdf[['Video #', 'Video Name', 'Cell', 'Track', 'Frame', 'x', 'y', 'Intensity',
                                    'isGap', 'GapFixed', 'Pass1', 'Pass2', 'Pass3', 'Bound']]
        output_tracks.append(trackdf)

    print_log('__________________________________________________')
    print_log('Complete: '
              '\n\t-> Frame Gap Filled:', counts_gap,
              '\n\t-> Events Separated:', counts_event,
              '\n\t-> Relabeling Performed:', counts_operations,
              '\n\t-> Fraction Filtered from total:', frames_filtered/frames_total,
              '\n\t-> Fraction Filtered from total:', frames_filtered_onlyBound / frames_total,
              )

    print_log('Saving to:', str(os.path.join(output_path, 'intermediates', 'gaps-and-fixes_decisions.csv')))
    pd.concat(output_tracks).to_csv(str(os.path.join(output_path, 'intermediates', 'gaps-and-fixes_decisions.csv')))

    return

'''
================================================================================================================
KEY FUNCTIONS
================================================================================================================
'''

# gap fixes, by frame
def process_gaps(track, pos, criteria):
    """
    Filling in the frame gaps based on the criteria provided.
    :param track: (list(list(float))) A track (list of spot information) separated and behavior classified.
    :param pos: (list(list(float))) position of spots by each frame, track with only position information.
    :param criteria: (lambda(int, int) -> int) function to determine gap behavior given the behaviors of neighboring frames.
    :return: count: number of gaps filled
    :return: result: processed track list
    :return: result_pos: processed pos list
    """
    # track[f][3] -> frame #
    # track[f][8] -> behavior
    result = []
    result_pos = []
    behaviors = [[0, 0, 0, 0], [1, 0, 1, 1], [1, 1, 0, 2]]
    count = 0
    f = 1
    result.append(np.array(list(track[0].copy()) + [0]))
    result_pos.append(pos[0])
    while f < len(track):
        if(track[f][3] - track[f-1][3] < 2):
            result.append(np.array(list(track[f].copy()) + [0]))
            result_pos.append(pos[f])
            f += 1
            continue
        fr = track[f-1][3] + 1
        bound = behaviors[criteria(track[f-1][8], track[f][8], track[f][3]-track[f-1][3] + 1)] # Behavior
        template = list(track[f-1][:3]).copy()

        # Gap filling
        while fr < track[f][3]:
            result.append(
                np.array(template.copy() + [fr, 0] + bound.copy() + [-1])
            )
            result_pos.append(np.array([-1, -1]))
            count += 1
            fr += 1
        result.append(np.array(list(track[f].copy()) + [0]))
        result_pos.append(pos[f])
        f += 1

    return count, result, result_pos

def pass_events(events, ori, sub, min_time):
    """
    relabel events with track length less than min time, then merge tracks with same behaviors
    :param events: (list(tuple(int, list(list(float)))) A record of every event in the form of (behavior, track).
    :param ori: (int) the behavior of interest of this pass
    :param sub: (lambda(int, int) -> int) The substituted behavior based on behaviors of neighboring frames.
    :param min_time: (int) minimum length in frames to consider a labelled event as invalid.
    :return: count: number of events relabelled.
    :return: result_events: all events processed.
    """
    count = 0
    behaviors = [[0, 0, 0, 0], [1, 0, 1, 1], [1, 1, 0, 2]]
    events = events.copy()
    for i in range(len(events)):
        bvr, event = events[i]
        event = event.copy()
        if not type(sub) == int:
            l = -1 if i-1 < 0 else events[i-1][0]
            r = -1 if i+1 >= len(events) else events[i+1][0]
        if not bvr == ori:
            continue
        if len(event) >= min_time:
            continue
        count += 1
        for f in range(len(event)):
            if type(sub) == int:
                event[f] = np.array(list(event[f][:5]) + behaviors[sub] + [event[f][9]])
            else:
                event[f] = np.array(list(event[f][:5]) + behaviors[sub(l, r)] + [event[f][9]])
        if type(sub) == int:
            events[i] = (sub, event)
        else:
            events[i] = (sub(l,r), event)
    result_events = []
    record_events = [events[0]]
    i = 1
    while i < len(events):
        if not events[i][0] == record_events[-1][0]:
            if len(record_events) == 1:
                result_events += record_events
                record_events = []
            else:
                result_events.append(merge_events(record_events))
                record_events = []
        record_events.append(events[i])
        i += 1
    if len(record_events) == 1:
        result_events += record_events
    else:
        result_events.append(merge_events(record_events))
    return count, result_events

'''
================================================================================================================
HELPER FUNCTIONS
================================================================================================================
'''

# Proportion calculation, weighted by behavior
def prop_counting(track, weights, countGaps=True):
    w1 = weights[0] # constricted
    w2 = weights[1] # strict
    gapCount = 0
    count1 = 0.0
    count2 = 0.0
    for spot in track:
        if not countGaps and spot[-1] == -1:
            gapCount += 1
            continue
        elif spot[-2] == 2:
            count2 += w2
        elif spot[-2] == 1:
            count1 += w1
    return (count1 + count2) / (len(track) - gapCount)


# Reconstruct Track from Events
def events_to_track(events):
    track = []
    for bvr, event in events:
        track += event
    return track


# merge events with the same behavior
def merge_events(events:list):
    event = []
    for i in range(len(events)):
        event += events[i][1]
    return (events[0][0], event)

# event separation, list(tuples(int, list(ndarray))))
def event_separation(track):
    # [f][8] -> bound
    count = 0
    events = []
    dtftransition = 0
    event = [track[0]]
    f = 1
    while f < len(track):

        if track[f][8] == 2 and (event[-1][8] == 1 or event[-1][8] == 0):
            dtftransition += 1

        if (track[f][8] == 1 or track[f][8] == 0) and event[-1][8] == 2:
            dtftransition += 1

        if track[f][8] == event[-1][8]:
            event.append(track[f])
        else:
            events.append((event[-1][8], event))
            count += 1
            event = [track[f]]
        f += 1
    if len(event) > 0:
        events.append((event[-1][8], event))
        count += 1
    return count, events, dtftransition

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

# write csv file at given path
def csv_write(path, data):
    with open(path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        for line in data:
            writer.writerow(line)
        file.close()

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
        prog='gaps-and-fixes',
        description='Further refining track behavior interpretation by interpolating the gap behavior.',
        epilog='Prior scripts: bound_classification.py')
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__ \n\t-> config_path: ' + str(args.config))
    main(args.config)
    print_log("--- %s seconds ---" % (time.time() - start_time))