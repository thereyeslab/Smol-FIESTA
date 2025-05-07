"""
Runnable Script if run as __main__
Calculate relevant track characteristics by tracing track behaviors and rebinding events.
- Determine rebinding events from the gap-fixed (or not) and classified track behaviors.
- Determine the difference in position of binding in rebinding events, classify them as same/different.
- Calculate frame-by-frame rebinding probability from *qualifying events, along with other parameters.
- Output rebinding events to Trackmate *.csv format for other tools.

Input:
    Either of:
        {csv_path}/{output_folder_name}/gaps-and-fixes_decisions.csv
        {csv_path}/{output_folder_name}/bound_decisions.csv
    Determined by parameter: {use_gap_fixed}

Output:
    {csv_path}/{output_folder_name}/RESULT_rebind.txt
    {csv_path}/{output_folder_name}/rebind-strict-event.csv
    {csv_path}/{output_folder_name}/rebind-strict-boundtime.csv
    {csv_path}/{output_folder_name}/rebind-flanked-strict-boundtime.csv
    {csv_path}/{output_folder_name}/rebind-AllDiffusion-time.csv
    {csv_path}/{output_folder_name}/rebind-strict-rebindingtime.csv
    {csv_path}/{output_folder_name}/SMAUG_REBINDING_SPOTS/strict_rebinds_spotsRebind.csv
    {csv_path}/{output_folder_name}/SMAUG_REBINDING_SPOTS/strict_rebinds_spotsSame.csv
    {csv_path}/{output_folder_name}/SMAUG_REBINDING_SPOTS/strict_rebinds_spotsDiff.csv
    {csv_path}/{output_folder_name}/SMAUG_REBINDING_SPOTS/strict_rebinds_spotsAll.csv

Parameters:
    rebind_distance_same: Determines rebinds to same particles if < parameter
    rebind_distance_diff: Determines rebinds to diff particles if > parameter
    min_time_bound_strict: min length (frames) of strict binding event.
    min_time_bound_constrained: min length (frames) of constrained diffusion event.
    min_time_rebinding_relaxed: min length (frames) of rebinding events considering relaxed criteria.
    min_time_rebinding_strict: min length (frames) of rebinding events considering strict criteria.
    min_time_diffusion: min length (frames) of free diffusion event
    min_time_diffusion_subsequent: minimum length (frames) required for subsequent event after diffusion.
    max_time_rebinding: max length (frames) after which rebinding will be counted as unsuccessful.
    max_time_constrained: max length (frames) after which rebinding will be counted as unsuccessful.

    conditional:
        use_gap_fixed: Use tracks processed by gaps_and_fixes.py, if False use tracks only processed by bound_classification.py instead.
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
    rebind_distance_same = configs['rebind-analysis'].get('rebind_distance_same', 'None')
    rebind_distance_diff = configs['rebind-analysis'].get('rebind_distance_diff', 'None')
    min_time_bound_strict = configs['rebind-analysis'].get('min_time_bound_strict', 'None')
    min_time_bound_constricted = configs['rebind-analysis'].get('min_time_bound_constrained', 'None')
    min_time_rebinding_relaxed = configs['rebind-analysis'].get('min_time_rebinding_relaxed', 'None')
    min_time_rebinding_strict = configs['rebind-analysis'].get('min_time_rebinding_strict', 'None')
    min_time_diffusion = configs['rebind-analysis'].get('min_time_diffusion', 'None')
    min_time_diffusion_subsequent = configs['rebind-analysis'].get('min_time_diffusion_subsequent', 'None')
    max_time_rebinding = configs['rebind-analysis'].get('max_time_rebinding', 'None')
    max_time_constrained = configs['rebind-analysis'].get('max_time_constrained', 'None')

    use_gap_fixed = configs['toggle']['use_gap_fixed']

    if isinstance(rebind_distance_same, str): rebind_distance_same = -1
    if isinstance(rebind_distance_diff, str): rebind_distance_diff = 100000
    if isinstance(min_time_bound_strict, str): min_time_bound_strict = 0
    if isinstance(min_time_bound_constricted, str): min_time_bound_constricted = 0
    if isinstance(min_time_rebinding_relaxed, str): min_time_rebinding_relaxed = 0
    if isinstance(min_time_rebinding_strict, str): min_time_rebinding_strict = 0
    if isinstance(min_time_diffusion, str): min_time_diffusion = 0
    if isinstance(min_time_diffusion_subsequent, str): min_time_diffusion_subsequent = 1
    if isinstance(max_time_rebinding, str): max_time_rebinding = 1000000
    if isinstance(max_time_constrained, str): max_time_constrained = 100000

    output_path = str(os.path.join(csv_path, configs['path']['output_folder_name']))
    log_file = str(os.path.join(output_path, 'LOG_rebind.txt'))
    log_result = str(os.path.join(output_path, 'RESULT_rebind.txt'))
    csv_result = str(os.path.join(output_path, 'RESULT_rebind.csv'))
    logging_setup(output_path, 'rebind')
    if not os.path.isdir(output_path):
        raise ValueError('Directory do not exist, please run track_sorting.py first.')

    tracks = pd.read_csv(str(os.path.join(output_path, ('gaps-and-fixes_decisions.csv' if use_gap_fixed else 'bound_decisions.csv'))))
    tracks = tracks.loc[:, ~tracks.columns.str.contains('^Unnamed')]

    headers = tracks[['Video #', 'Cell', 'Track']].to_numpy()
    tracks = slice_tracks(tracks, headers)
    headers = np.unique(headers, axis=0)

    rebind_relaxed = []
    rebind_relaxed_spots_same = []
    rebind_relaxed_spots_diff = []
    rebind_relaxed_spots_all = []
    rebind_relaxed_time_all = []
    rebind_relaxed_unsuccessful = 0
    bound_constricted = []
    bound_constricted_record = []

    bound_flanked_relaxed_record = []

    rebind_relaxed_spots_entiretrack = []

    rebind_strict = []
    rebind_strict_spots_same = []
    rebind_strict_spots_diff = []
    rebind_strict_spots_all = []
    rebind_strict_time_all = []
    rebind_strict_unsuccessful = 0
    bound_strict = []
    bound_strict_record = []

    bound_flanked_strict_record = []

    rebind_strict_spots_entiretrack = []

    constrained_dest = np.array([0, 0])

    fast_diffusion_dest = np.array([0, 0])
    fast_diffusion_dest_strict = np.array([0, 0])
    fast_diffusion_time = []

    all_diffusion_dest = np.array([0, 0])
    all_diffusion_dest_strict = np.array([0, 0])
    all_diffusion_time = []
    diftime_record = []

    proportion_count = np.array([0, 0, 0])

    abound_dest = np.array([0, 0])

    for i in range(len(tracks)):
        header = headers[i]
        track = list(tracks[i][['Frame', 'x', 'y', 'Bound']].to_numpy())

        # Relaxed
        rb, rb_us, rb_same, rb_diff, rb_all, time_all = (
            rebind_record_proximity(track, rebind_distance_same, rebind_distance_diff, lambda x: not x < 1, min_time_rebinding_relaxed, max_time_rebinding))
        bd = bound_record(track, lambda x: x == 1, min_time_bound_constricted)
        if(len(rb) > 0):
            for j in range(len(rb)):
                rb[j] = list(header) + rb[j]
            rebind_relaxed += rb
            rebind_relaxed_spots_entiretrack.append([track.copy()])
        rebind_relaxed_unsuccessful += rb_us
        if(len(rb_same) > 0):
            rebind_relaxed_spots_same.append(rb_same)
        if(len(rb_diff) > 0):
            rebind_relaxed_spots_diff.append(rb_diff)
        if(len(rb_all) > 0):
            rebind_relaxed_spots_all.append(rb_all)
        if(len(bd) > 0):
            bound_constricted += bd
        j = 1
        for bdframe in bd:
            bound_constricted_record.append(list(header.copy()) + [j, bdframe])
            j += 1
        j = 1
        for rbtime in time_all:
            rebind_relaxed_time_all.append(list(header.copy()) + [j, rbtime])
            j += 1

        # Bound time calculations, now must be flanked by diffusions
        _, __, ___, ____, _____, bd1 = (
            rebind_record_proximity(track, rebind_distance_same, rebind_distance_diff, lambda x: not x > 0,
                                    min_time_bound_constricted, max_time_rebinding))
        j = 1
        for bdtime in bd1:
            bound_flanked_relaxed_record.append(list(header.copy()) + [j, bdtime])
            j += 1

        # Strict
        rb, rb_us, rb_same, rb_diff, rb_all, time_all = (
            rebind_record_proximity(track, rebind_distance_same, rebind_distance_diff, lambda x: not x < 2, min_time_rebinding_strict, max_time_rebinding))
        bd = bound_record(track, lambda x: x == 2, min_time_bound_strict)
        if(len(rb) > 0):
            for j in range(len(rb)):
                rb[j] = list(header) + rb[j]
            rebind_strict += rb
            rebind_strict_spots_entiretrack.append([track.copy()])
        rebind_strict_unsuccessful += rb_us
        if(len(rb_same) > 0):
            rebind_strict_spots_same.append(rb_same)
        if(len(rb_diff) > 0):
            rebind_strict_spots_diff.append(rb_diff)
        if(len(rb_all) > 0):
            rebind_strict_spots_all.append(rb_all)
        if(len(bd) > 0):
            bound_strict += bd
        j = 1
        for bdframe in bd:
            bound_strict_record.append(list(header.copy()) + [j, bdframe])
            j += 1
        j = 1
        for rbtime in time_all:
            rebind_strict_time_all.append(list(header.copy()) + [j, rbtime])
            j += 1

        ########

        bf_time, bf_counts = diffusion_record(track, lambda x: x < 2, min_time_bound_strict,
                                              min_time_diffusion_subsequent)
        abound_dest = np.add(abound_dest, bf_counts)
        ########

        # Bound time calculations, now must be flanked by diffusions
        _, __, ___, ____, _____, bd1 = (
            rebind_record_proximity(track, rebind_distance_same, rebind_distance_diff, lambda x: not x > 1,
                                    min_time_bound_strict, max_time_rebinding))
        j = 1
        for bdtime in bd1:
            bound_flanked_strict_record.append(list(header.copy()) + [j, bdtime])
            j += 1

        # constrained diffusion
        constrained_dest = np.add(constrained_dest, constrained_record(track, min_time_bound_constricted, min_time_bound_strict, max_time_constrained))

        # fast diffusion
        df_time, df_counts = diffusion_record(track, lambda x: x > 0, min_time_diffusion, min_time_diffusion_subsequent)
        fast_diffusion_time += df_time
        fast_diffusion_dest = np.add(fast_diffusion_dest, df_counts)
        if len(bd) > 0:
            fast_diffusion_dest_strict = np.add(fast_diffusion_dest_strict, df_counts)

        # all diffusion
        df_time, df_counts = diffusion_record(track, lambda x: x > 1, min_time_diffusion,
                                              min_time_diffusion_subsequent)
        all_diffusion_time += df_time
        all_diffusion_dest = np.add(all_diffusion_dest, df_counts)
        if len(bd) > 0:
            all_diffusion_dest_strict = np.add(all_diffusion_dest_strict, df_counts)
        j = 1
        for dfframe in df_time:
            diftime_record.append(list(header.copy()) + [j, dfframe])
            j += 1

        # proportion count
        p_counts = label_count(track)
        proportion_count = np.add(proportion_count, p_counts)

    # for a csv output of the results
    output_result = []
    print_log('[Analysis]')

    # binding events output
    print_log('__________Bound__________')
    output_result.append('Bound')

    print_log('Constrained Diffusion Time (Frame):')
    print_log('->', str(pd.Series(bound_constricted).describe()).replace('\n','\n-> '),'\n')
    output_result.append(
        'Constrained Diffusion Time (Frame),' + ','.join(str(pd.Series(bound_constricted).describe()).split()[:-2]))
    print_log('Strict Bound Time (Frame):')
    print_log('->', str(pd.Series(bound_strict).describe()).replace('\n','\n-> '), '\n')
    output_result.append(
        'Strict Bound Time (Frame),' + ','.join(str(pd.Series(bound_strict).describe()).split()[:-2]))

    print_log('\n______Bound_Flanked_______')
    print_log('# Note: This is bound time flanked by diffusion.')
    print_log('\nStrict Bound time:')
    print_log('->', str(pd.Series([x[4] for x in bound_flanked_strict_record]).describe()).replace('\n', '\n-> '), '\n')
    output_result.append(
        'Flanked Bound Time (Frame),' + ','.join(str(pd.Series([x[4] for x in bound_flanked_strict_record]).describe()).split()[:-2]))

    output_result.append('')

    print_log('\nBound to diffusion by Frame')
    print_log('-> Count of Bound to Bound:', abound_dest[0])
    print_log('-> Count of Bound to Diffusion:', abound_dest[1])

    if abound_dest[1] == 0 or abound_dest[0] == 0:
        op = 0
    else:
        op = float(abound_dest[1]) / float(abound_dest[0] + abound_dest[1])

    print_log('-> Probability of Bound to diffusion:', op)
    output_result.append(
        'Bound Transitions (by Frame),B->B,' + str(abound_dest[0]) + ',B->D,' + str(abound_dest[1]) + ',P(B->D),' + str(op))

    output_result.append('')

    # rebinding events output
    print_log('__________Rebind_________')
    output_result.append('Rebind')

    print_log('\n Strict bindings Rebind Time (Frame):')
    print_log('->', str(pd.Series([x[4] for x in rebind_strict_time_all]).describe()).replace('\n', '\n-> '), '\n')
    output_result.append(
        'Strict Bindings Rebind Time (Frame),' + ','.join(str(pd.Series([x[4] for x in rebind_strict_time_all]).describe()).split()[:-2]))

    print_log('Strict bindings Rebind Probability:')
    print_log('-> Successful:', len(rebind_strict))
    print_log('-> Unsuccessful:', rebind_strict_unsuccessful)

    if len(rebind_strict) == 0 or rebind_strict_unsuccessful == 0:
        op1 = 0
    else:
        op1 = float(len(rebind_strict)) / float(len(rebind_strict) + rebind_strict_unsuccessful)

    print_log('-> Probability', op1, '\n' )
    output_result.append(
        'Rebinding (Events),Success,' + str(len(rebind_strict)) + ',Fail,' +
        str(rebind_strict_unsuccessful) + ',P(Success),' + str(op1))

    output_result.append('')

    # constrained diffusion output
    print_log('_______Constrained________')
    output_result.append('Constrained')

    print_log('Count of Constrained to Diffusion:', constrained_dest[0])
    print_log('Count of Constrained to Bound:', constrained_dest[1])
    print_log('Probability of Constrained to Bound:', 0.000001+float(constrained_dest[1]) / 0.000001+(float(constrained_dest[0] + constrained_dest[1])))
    print_log('')
    output_result.append(
        'Constrained Diffusion Transitions (by Events),->Diffusion,' + str(constrained_dest[0]) + ',->Bound,' +
        str(constrained_dest[1]) + ',P(->Bound),' + str(0.0000001+(float(constrained_dest[1]))/0.000001 +(float(constrained_dest[0] + constrained_dest[1])+0.0000001)))

    output_result.append('')

    # fast diffusion output
    print_log('______Fast_Diffusion______')
    output_result.append('Fast Diffusion')
    print_log('\nFast Diffusion average time (Frame):')
    print_log('->', str(pd.Series(fast_diffusion_time).describe()).replace('\n', '\n-> '), '\n')
    output_result.append(
        'Fast Diffusion Time (Frame),' + ','.join(str(pd.Series(fast_diffusion_time).describe()).split()[:-2]))

    output_result.append('')

    # All diffusion events output
    print_log('______All_Diffusion_______')
    output_result.append('All Diffusion (Fast and Constrained)')
    print_log('All Diffusion by Frame')
    print_log('-> Count of Diffusion to Diffusion:', all_diffusion_dest[0])
    print_log('-> Count of Diffusion to Bound:', all_diffusion_dest[1])
    print_log('-> Probability of Diffusion to Bound:',
              float(all_diffusion_dest[1]) / (float(all_diffusion_dest[0] + all_diffusion_dest[1]))+.000000001)
    output_result.append(
        'Diffusion(All) Transitions (by Frame),->Diffusion,' + str(all_diffusion_dest[0]) + ',->Bound,' +
        str(all_diffusion_dest[1]) + ',P(->Bound),' +
        str(float(all_diffusion_dest[1]) / float(all_diffusion_dest[0] + all_diffusion_dest[1])))

    print_log('\nAll Diffusion (with strict binding event in track) by Frame')
    print_log('-> Count of Diffusion to Diffusion:', all_diffusion_dest_strict[0])
    print_log('-> Count of Diffusion to Bound:', all_diffusion_dest_strict[1])

    if all_diffusion_dest_strict[0] == 0 or all_diffusion_dest_strict[1] == 0:
        op2 = 0
    else:
        op2 = float(all_diffusion_dest_strict[1]) / (float(all_diffusion_dest_strict[0] + all_diffusion_dest_strict[1])+0.000000001)

    print_log('-> Probability of Diffusion to Bound:', op2)
    output_result.append(
        'Diffusion(with strict binding event in track) Transitions (by Frame),->Diffusion,' +
        str(all_diffusion_dest_strict[0]) + ',->Bound,' +
        str(all_diffusion_dest_strict[1]) + ',P(->Bound),' +
        str(op2))

    print_log('\nAll Diffusion average time (Frame):')
    print_log('->', str(pd.Series(all_diffusion_time).describe()).replace('\n', '\n-> '), '\n')
    output_result.append(
        'Diffusion(All) Time (Frame),' + ','.join(str(pd.Series(all_diffusion_time).describe()).split()[:-2]))

    output_result.append('')

    # All events count by frame output
    print_log('____Counted_All_Events____')
    output_result.append('All Events')

    print_log('Count of all frames:', np.sum(proportion_count))
    output_result.append('All (Frame),count,' + str(np.sum(proportion_count)))

    print_log('Count of fast diffusion: ', proportion_count[0], '-> Proportion:',
              float(proportion_count[0]) / np.sum(proportion_count))
    output_result.append('Fast Diffusion (Frame),count,' + str(proportion_count[0]) +
                         ',proportion,' + str(float(proportion_count[0]) / np.sum(proportion_count)))

    print_log('Count of constrained diffusion: ', proportion_count[1], '-> Proportion:',
              float(proportion_count[1]) / np.sum(proportion_count))
    output_result.append('Constrained Diffusion (Frame),count,' + str(proportion_count[1]) +
                         ',proportion,' + str(float(proportion_count[1]) / np.sum(proportion_count)))

    if proportion_count[2] == 0:
        op3 = 0
    else:
        op3 = float((proportion_count[2]) / np.sum(proportion_count))

    print_log('Count of strict binding: ', proportion_count[2], '-> Proportion:', op3)
    output_result.append('Bound (Frame),count,' + str(proportion_count[2]) +
                         ',proportion,' + str(op3))

    # output, truncate log_RESULT
    with open(log_file) as fin, open(log_result, 'w') as fout:
        active = False
        for line in fin:
            if '[Analysis]' in line:
                active = True
            if active:
                fout.write(line)

    # output, csv_RESULT
    with open(csv_result, "w") as csv_file:
        for line in output_result:
            csv_file.write(line)
            csv_file.write('\n')

    # outputs
    rebind_strict_spots_all = event_format_trackmate(rebind_strict_spots_all)
    rebind_strict_spots_same = event_format_trackmate(rebind_strict_spots_same)
    rebind_strict_spots_diff = event_format_trackmate(rebind_strict_spots_diff)
    rebind_strict_spots_entiretrack = event_format_trackmate(rebind_strict_spots_entiretrack)

    smaug_path = str(os.path.join(output_path, 'SMAUG_REBINDING_SPOTS'))
    try:
        shutil.rmtree(smaug_path)
        os.mkdir(smaug_path)
    except:
        os.mkdir(smaug_path)
    csv_write(str(os.path.join(smaug_path, 'strict_rebinds_spotsRebind.csv')), rebind_strict_spots_all)
    csv_write(str(os.path.join(smaug_path, 'strict_rebinds_spotsSame.csv')), rebind_strict_spots_same)
    csv_write(str(os.path.join(smaug_path, 'strict_rebinds_spotsDiff.csv')), rebind_strict_spots_diff)
    csv_write(str(os.path.join(smaug_path, 'strict_rebinds_spotsAll.csv')), rebind_strict_spots_entiretrack)

    rebind_columns = ['Video #', 'Cell', 'Track', 'From', 'To', 'Time', 'Speed', 'Distance', 'x1', 'y1', 'x2', 'y2']
    rebind_strict = pd.DataFrame(rebind_strict, columns=rebind_columns).astype({'Time': 'int'})
    rebind_strict.to_csv(str(os.path.join(output_path, 'rebind-strict-event.csv')))

    boundtime_columns = ['Video #', 'Cell', 'Track', 'Event', 'Bound Time']
    bound_constricted_record = pd.DataFrame(bound_constricted_record, columns=boundtime_columns).astype({'Bound Time': 'int'})
    bound_strict_record = pd.DataFrame(bound_strict_record, columns=boundtime_columns).astype({'Bound Time': 'int'})
    bound_constricted_record.to_csv(str(os.path.join(output_path, 'rebind-constrained-DiffusionTime.csv')))
    bound_strict_record.to_csv(str(os.path.join(output_path, 'rebind-strict-boundtime.csv')))

    bound_flanked_strict_record = pd.DataFrame(bound_flanked_strict_record, columns=boundtime_columns).astype({'Bound Time': 'int'})
    bound_flanked_strict_record.to_csv(str(os.path.join(output_path, 'rebind-flanked-strict-boundtime.csv')))

    diftime_columns = ['Video #', 'Cell', 'Track', 'Event', 'Diffusion Time']
    diftime_record = pd.DataFrame(diftime_record, columns=diftime_columns).astype(
        {'Diffusion Time': 'int'})
    diftime_record.to_csv(str(os.path.join(output_path, 'rebind-AllDiffusion-time.csv')))

    rbtime_columns = ['Video #', 'Cell', 'Track', 'Event', 'Rebinding Time']
    rbtime_strict_record = pd.DataFrame(rebind_strict_time_all, columns=rbtime_columns).astype(
        {'Rebinding Time': 'int'})
    rbtime_strict_record.to_csv(str(os.path.join(output_path, 'rebind-strict-rebindingtime.csv')))
    # Transition Matrix calculation
    transition_matrices = calculate_transition_matrices(tracks)

    # Save corrected CSV output
    transition_csv = os.path.join(output_path, 'transition_matrices_counts_and_proportions.csv')
    transition_matrices.to_csv(transition_csv)

    # Log updated matrices
    print_log('Corrected transition matrices (frame-by-frame counts & proportions):')
    print_log(f'\n{transition_matrices}')

    return

'''
================================================================================================================
KEY FUNCTION
================================================================================================================
'''

def rebind_record_proximity(track, rebind_distance_same, rebind_distance_diff, criteria, min_time, max_time):
    """
    Process and trace tracks for detecting rebinding behaviors of particle spots. Also tabulates relevant information.
    :param track: (list(list(float))) Track, list of spot information, labelled by the behaviors for each spot.
    :param rebind_distance_same: (int) distance in pixels for the beginning and end of rebinding events to be considered same positions.
    :param rebind_distance_diff: (int) distance in pixels for the beginning and end of rebinding events to be considered different positions.
    :param criteria: (lambda(int) -> bool) criteria for binding behavior, * >= 2 for strict binding.
    :param min_time: (int) minimum rebinding event duration.
    :param max_time: (int) maximum rebinding event duration, after which the event is considered to be failed.
    :return: rebinds: (list(list(float))) all rebinding events, with added trace information.
    :return: unsuccessful: (int) number of unsuccessful rebinding events.
    :return: events_same: (list(list(float))) rebinding events start and end in the same position.
    :return: events_diff: (list(list(float))) rebinding events start and end in different positions.
    :return: events_all: (list(list(float))) all rebinding events, same format as the input tracks.
    :return: time_all: (list(int)) rebinding time in frames for all rebinding events.
    """
    rebinds = []
    event = []
    events_same = []
    events_diff = []
    events_all = []
    time_all = []
    active = False
    record_pos = [track[0][1], track[0][2]]
    record_f = 0
    unsuccessful = 0
    f = 0

    while (f < len(track)):
        if (len(event) > 0 and criteria(track[f][3])):
            pos = [track[f][1], track[f][2]]
            dist = distance(pos, record_pos)
            table = rebind_tabulate(event.copy(), 0, 0) # just to get the time
            if(table[2] < min_time): # min_time threshold
                event = []
                time_all.append(table[2])
            elif(table[2] > max_time):
                event = []
                time_all.append(table[2])
                unsuccessful += 1
            else:
                time_all.append(table[2])
                if (dist >= rebind_distance_diff):
                    prev, nxt = 1, 2
                    events_diff.append(event.copy())
                elif (dist <= rebind_distance_same):
                    prev, nxt = 1, 1
                    events_same.append(event.copy())
                else:
                    prev, nxt = 1, 2
                events_all.append(event.copy())
                time_int = rebind_trace_avg(track, f, criteria, 1)[0]
                rebinds.append(
                    rebind_tabulate(event.copy(), prev, nxt) + [dist] +
                    rebind_trace_avg(track, record_f - 1, criteria, -1)[1] +
                    rebind_trace_avg(track, f, criteria, 1)[1]
                )
                event = []
        if criteria(track[f][3]):
            active = True
            record_pos = [track[f][1], track[f][2]]
        elif (active):
            if len(event) == 0:
                record_f = f
            event.append(track[f])
        f += 1

    # unsuccessful event
    unsuccessful += 1 if (len(event) > 0) else 0
    return rebinds, unsuccessful, events_same, events_diff, events_all, time_all

'''
================================================================================================================
BEHAVIOR PROCESSING
================================================================================================================
'''

# Process diffusion behaviors in tracks by tracing, records diffusion time and endpoint behavior (sufficently bound or not)
def diffusion_record(track, criteria, min_time, min_time_bound):
    counts = [0, 0]
    event = []
    event2 = []
    diffusion_time = []
    f = 0
    while f < len(track):
        if not criteria(track[f][3]):
            event.append(track[f])
        else:
            if len(event) > 0:
                time_int = 1 if len(event) == 1 else event[-1][0] - event[0][0] + 1
                if time_int < min_time:
                    event = []
                    continue
                diffusion_time.append(time_int)
                counts[0] += len(event) - 1
                record = track[f][3]
                i = f
                while i < len(track) and criteria(track[i][3]):
                    event2.append(track[i])
                    i += 1
                time_int2 = 1 if len(event2) == 1 else event2[-1][0] - event2[0][0] + 1
                if time_int2 < min_time_bound:
                    counts[0] += 1
                else:
                    counts[1] += 1
                event2 = []
                event = []
        f += 1
    if len(event) > 0:
        time_int = 1 if len(event) == 1 else event[-1][0] - event[0][0] + 1
        if time_int >= min_time:
            diffusion_time.append(time_int)
            counts[0] += len(event) - 1
    return diffusion_time, np.array(counts)

# Process constrained-diffusion behaviors in tracks by tracing, records counts of endpoint behavior (whether diffusion or bound)
def constrained_record(track, min_time_constrained, min_time_bound, max_time_constrained):
    counts = [0, 0]
    f = 0
    event = []
    event2 = []
    while f < len(track):
        if track[f][3] == 1:
            event.append(track[f])
        else:
            if(len(event) > 0):
                time_int = 1 if len(event) == 1 else event[-1][0] - event[0][0] + 1
                record = track[f][3]
                i = f
                while i < len(track) and track[i][3] == record:
                    event2.append(track[i])
                    i += 1
                time_int2 = 1 if len(event2) == 1 else event2[-1][0] - event2[0][0] + 1
                if(time_int >= min_time_constrained and time_int <= max_time_constrained and
                        time_int2 >= min_time_bound):
                    counts[0 if track[f][3] == 0 else 1] += 1
                event2 = []
                event = []
        f += 1

    return np.array(counts)

# Process bound behaviors in tracks by tracking, records time for each bound event.
def bound_record(track, criteria, min_time):
    bound = []
    record = track[0][3]
    event = []
    f = 0
    while f < len(track):
        if record == track[f][3]:
            if criteria(track[f][3]):
                event.append(track[f])
        else:
            if(len(event) > 0):
                time = 1 if len(event) == 1 else int(event[-1][0] - event[0][0] + 1)
                if(time < min_time):
                    event = []
                else:
                    bound.append(event.copy())
                    event = []
            if(criteria(track[f][3])):
                event.append(track[f])
            record = track[f][3]
        f += 1
    if(len(event) > 0):
        time = 1 if len(event) == 1 else int(event[-1][0] - event[0][0] + 1)
        if (time < min_time):
            event = []
        else: bound.append(event.copy())

    result = []
    for event in bound:
        if(len(event) == 1):
            result.append(1)
        else:
            result.append(int(event[-1][0] - event[0][0] + 1))
    return result


'''
================================================================================================================
HELPER FUNCTIONS
================================================================================================================
'''

# returns counts for each track label (diffusion, constrained, bound)
def label_count(track):
    result = np.array([0, 0, 0])
    labels = np.array(track)[:, 3]
    unique, counts = np.unique(labels, return_counts=True)
    unique = unique.astype('int')
    for i in range(unique.shape[0]):
        result[unique[i]] = counts[i]
    return result

# format each event to Trackmate format, but loses intensity information
def event_format_trackmate(events):
    formatted = []
    i = 1
    for track in events:
        for event in track:
            for spot in event:
                if spot[1] < 0 or spot[2] < 0:
                    continue
                formatted.append([i, spot[0], spot[1], spot[2], 10000])
            i += 1
    return formatted

# Find the average position of a track
def rebind_trace_avg(track, sframe, criteria, dir):
    f = sframe
    x = []
    y = []
    event = []
    while f >= 0 and f < len(track) and criteria(track[f][3]):
        event.append(track[f])
        x.append(track[f][1])
        y.append(track[f][2])
        f += dir
    time_int = 1 if len(event) == 1 else int(event[-1][0] - event[0][0] + 1)
    return time_int, [np.mean(x), np.mean(y)]

# Iterate through a track segment and calculate distance travelled by spot.
def rebind_tabulate(segment, prev, nxt):
    frames = [s[0] for s in segment]
    rebinding_time = max(frames) - min(frames) + 1
    distances = [0]
    if (len(segment) > 1):
        for i in range(1, len(segment)):
            distances.append(distance([segment[i - 1][1], segment[i - 1][2]], [segment[i][1], segment[i][2]]))
    return [prev, nxt, rebinding_time, sum(distances) / rebinding_time if (len(segment) > 1) else 1]

def distance(p1, p2):
    return np.sqrt(np.power(p1[0]-p2[0], 2) + np.power(p1[1] - p2[1], 2))

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

def csv_write(path, data):
    with open(path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        for line in data:
            writer.writerow(line)
        file.close()

def calculate_transition_matrices(tracks):
    states = {0: 'F.dif', 1: 'C.Dif', 2: 'Bound'}

    # Absolute transition count matrix
    abs_matrix = pd.DataFrame(
        np.zeros((3, 3), dtype=int),
        index=[f'From_{states[i]}' for i in range(3)],
        columns=[f'To_{states[i]}' for i in range(3)]
    )

    # Count transitions frame-by-frame
    for track in tracks:
        behavior_sequence = track['Bound'].values
        frame_numbers = track['Frame'].values

        for i in range(1, len(behavior_sequence)):
            prev_frame = frame_numbers[i-1]
            current_frame = frame_numbers[i]

            # Check consecutive frames
            if current_frame == prev_frame + 1:
                from_state = behavior_sequence[i-1]
                to_state = behavior_sequence[i]
                abs_matrix.iloc[from_state, to_state] += 1

    # Calculate proportions row-wise
    proportion_matrix = abs_matrix.div(abs_matrix.sum(axis=1), axis=0).fillna(0)

    # Combine matrices side-by-side
    combined_matrix = pd.concat(
        [abs_matrix, proportion_matrix],
        axis=1,
        keys=['Absolute Counts', 'Proportions']
    )

    return combined_matrix



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