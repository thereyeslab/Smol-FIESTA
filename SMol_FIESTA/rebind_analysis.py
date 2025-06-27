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
    {csv_path}/{output_folder_name}/rebind-AllDiffusion-time.csv
    {csv_path}/{output_folder_name}/rebind-strict-rebindingtime.csv

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
    use_gap_fixed = configs['toggle']['use_gap_fixed']
    # config defaults
    rebind_defaults = {
        'rebind_distance_same': -1,
        'rebind_distance_diff': 100000,
        'min_time_bound_strict': 0,
        'min_time_bound_constrained': 0,
        'min_time_rebinding_relaxed': 0,
        'min_time_rebinding_strict': 0,
        'min_time_diffusion': 0,
        'min_time_diffusion_subsequent': 1,
        'max_time_rebinding': 1000000,
        'max_time_constrained': 100000
    }


    rebind_cfg = {**rebind_defaults, **configs.get('rebind-analysis', {})}

    # unpack configs
    rebind_distance_same = rebind_cfg['rebind_distance_same']
    rebind_distance_diff = rebind_cfg['rebind_distance_diff']
    min_time_bound_strict = rebind_cfg['min_time_bound_strict']
    min_time_bound_constricted = rebind_cfg['min_time_bound_constrained']
    min_time_rebinding_relaxed = rebind_cfg['min_time_rebinding_relaxed']
    min_time_rebinding_strict = rebind_cfg['min_time_rebinding_strict']
    min_time_diffusion = rebind_cfg['min_time_diffusion']
    min_time_diffusion_subsequent = rebind_cfg['min_time_diffusion_subsequent']
    max_time_rebinding = rebind_cfg['max_time_rebinding']
    max_time_constrained = rebind_cfg['max_time_constrained']

    output_path = str(os.path.join(csv_path, configs['path']['output_folder_name']))
    log_file = str(os.path.join(output_path, 'logs', 'LOG_rebind.txt'))
    log_result = str(os.path.join(output_path, 'RESULT_rebind.txt'))
    csv_result = str(os.path.join(output_path, 'RESULT_rebind.csv'))
    logging_setup(output_path, 'rebind')
    if not os.path.isdir(output_path):
        raise ValueError('Directory do not exist, please run track_sorting.py first.')

    tracks = pd.read_csv(str(os.path.join(output_path, 'intermediates', ('gaps-and-fixes_decisions.csv' if use_gap_fixed else 'bound_decisions.csv'))))
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


    rebind_relaxed_spots_entiretrack = []

    rebind_strict = []
    rebind_strict_spots_same = []
    rebind_strict_spots_diff = []
    rebind_strict_spots_all = []
    rebind_strict_time_all = []
    rebind_strict_unsuccessful = 0
    bound_strict = []
    bound_strict_record = []
    fastdiff_headers = []

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
        segments_bound_strict      = extract_segments(track, target_label=2)
        segments_bound_constricted = extract_segments(track, target_label=1)
        segments_fast_diffusion    = extract_segments(track, target_label=0)
        segments_all_diffusion     = extract_segments(track, target_label=3)  # or adjust as needed
        segments_rebind_strict     = extract_segments_from_rebind(track, time_all)


        # 2) Constrained-Diffusion events
        j = 1
        for segment in segments_bound_constricted:
            start_frame = segment[0][0]
            end_frame   = start_frame + len(segment) - 1
            bound_constricted_record.append(list(header.copy()) + [j, start_frame, end_frame])
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
        # 1) Strict-Bound events
        j = 1
        for segment in segments_bound_strict:
            start_frame = segment[0][0]
            end_frame   = start_frame + len(segment) - 1
            bound_strict_record.append(list(header.copy()) + [j, start_frame, end_frame])
            j += 1

        # 3) Fast-Diffusion events
        j = 1
        for segment in segments_fast_diffusion:
            start_frame = segment[0][0]
            end_frame   = start_frame + len(segment) - 1
            fastdiff_headers.append(list(header.copy()) + [j, start_frame, end_frame])
            j += 1

        ########

        bf_time, bf_counts = diffusion_record(track, lambda x: x < 2, min_time_bound_strict,
                                              min_time_diffusion_subsequent)
        abound_dest = np.add(abound_dest, bf_counts)
        ########


        # constrained diffusion
        constrained_dest = np.add(constrained_dest, constrained_record(track, min_time_bound_constricted, min_time_bound_strict, max_time_constrained))

        # fast diffusion
        df_time, df_counts = diffusion_record(track, lambda x: x > 0, min_time_diffusion, min_time_diffusion_subsequent)
        fast_diffusion_time += df_time

        for ft in df_time:
            fastdiff_headers.append(list(header.copy()) + [j, ft])
            j += 1

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
        segments_all_diffusion = extract_segments(track, exclude_labels={1, 2})
        # Rebinding: diffusion segments flanked by strict-bound (2)
        j_rb = 1
        f = 1
        while f < len(track) - 1:
            if track[f - 1][3] == 2 and track[f][3] in [0, 1]:
                start_frame = track[f][0]
                while f < len(track) and track[f][3] in [0, 1]:
                    f += 1
                if f < len(track) and track[f][3] == 2:
                    end_frame = track[f - 1][0]  # last diffusion frame
                    rebind_strict_time_all.append(list(header.copy()) + [j_rb, start_frame, end_frame])
                    j_rb += 1
            else:
                f += 1

        # SearchTime extraction (label 0 or 1, flanked by 2 or edges)

        j_st = 1
        in_segment = False
        segment_start = None

        for idx, spot in enumerate(track):
            label = int(spot[3])
            frame = int(spot[0])

            if label in [0, 1]:
                if not in_segment:
                    in_segment = True
                    segment_start = frame
            elif label == 2 and in_segment:
                # We were in a diffusion segment, now it ends
                diftime_record.append(list(header.copy()) + [j_st, segment_start, track[idx - 1][0]])
                j_st += 1
                in_segment = False

        # If the track ends while still in a diffusion segment
        if in_segment:
            diftime_record.append(list(header.copy()) + [j_st, segment_start, track[-1][0]])

        # proportion count
        p_counts = label_count(track)
        proportion_count = np.add(proportion_count, p_counts)

    # for a csv output of the results
    output_result = []
    print_log('[Analysis]')

    # binding events output
    print_log('_______________________________Event Mean Time_______________________________')
    stats_strict      = pd.Series(bound_strict).describe()
    stats_constrained = pd.Series(bound_constricted).describe()
    stats_fast        = pd.Series(fast_diffusion_time).describe()
    stats_search      = pd.Series(all_diffusion_time).describe()
    stats_rebind      = pd.Series([x[4] for x in rebind_strict_time_all]).describe()

    summary_table = pd.DataFrame({
        'Bound':      stats_strict,
        'C. Diffusion':  stats_constrained,
        'Fast Diffusion':         stats_fast,
        'Diffusion Time':            stats_search,
        'Rebind Time':            stats_rebind
    })

    print_log(summary_table, '\n')

    # Determine which legend to use based on max_time_rebinding
    if max_time_rebinding < 1000000:
        print_log(f"Rebinding probability after {max_time_rebinding} frames:")
    else:
        print_log("Probability of Successful Binding:")

    print_log("-> Successful:", len(rebind_strict))
    print_log("-> Unsuccessful:", rebind_strict_unsuccessful)

    if len(rebind_strict) == 0 or rebind_strict_unsuccessful == 0:
        op1 = 0
    else:
        op1 = float(len(rebind_strict)) / (len(rebind_strict) + rebind_strict_unsuccessful)

    print_log("-> Probability", op1, "\n")

    output_result.append(
        'Rebinding (Events),Success,' + str(len(rebind_strict)) + ',Fail,' +
        str(rebind_strict_unsuccessful) + ',P(Success),' + str(op1))

    output_result.append('')

    # All events count by frame output (handles zero‐event cases)
    total_frames = np.sum(proportion_count)

    # Safely extract each count (in case the list is shorter than expected)
    fast_count      = proportion_count[0] if len(proportion_count) > 0 else 0
    constrained_count = proportion_count[1] if len(proportion_count) > 1 else 0
    bound_count     = proportion_count[2] if len(proportion_count) > 2 else 0

    # Compute proportions only if total_frames > 0, else default to 0.0
    if total_frames > 0:
        prop_fast       = float(fast_count) / total_frames
        prop_constrained = float(constrained_count) / total_frames
        prop_bound      = float(bound_count) / total_frames
    else:
        prop_fast = prop_constrained = prop_bound = 0.0

    # Build a DataFrame for printing
    df_counts = pd.DataFrame({
        'Event': [ 'F. Diffusion', 'C. Diffusion', 'Bound'],
        'Count': [ fast_count, constrained_count, bound_count],
        'Proportion': [ prop_fast, prop_constrained, prop_bound]
    })

    print_log('_______Event proportions_______')
    print_log(df_counts.to_string(index=False), '\n')

    # If you append to output_result, use the same safe proportions:
    output_result.append('All Events')
    output_result.append(f'Fast Diffusion,count,{fast_count},proportion,{prop_fast}')
    output_result.append(f'C. Diffusion,count,{constrained_count},proportion,{prop_constrained}')
    output_result.append(f'Strict Binding,count,{bound_count},proportion,{prop_bound}')

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

    rebind_columns = ['Video #', 'Cell', 'Track', 'From', 'To', 'Time', 'Speed', 'Distance', 'x1', 'y1', 'x2', 'y2']
    rebind_strict = pd.DataFrame(rebind_strict, columns=rebind_columns).astype({'Time': 'int'})
    rebind_strict.to_csv(str(os.path.join(output_path, 'rebind-strict-event.csv')))

    boundtime_columns = ['Video #', 'Cell', 'Track', 'Event', 'StartFrame', 'EndFrame']
    diftime_columns = ['Video #', 'Cell', 'Track', 'Event', 'StartFrame', 'EndFrame']
    fastdiff_columns = ['Video #', 'Cell', 'Track', 'Event', 'StartFrame', 'EndFrame']
    rbtime_columns = ['Video #', 'Cell', 'Track', 'Event', 'StartFrame', 'EndFrame']

    df_constrained = pd.DataFrame(bound_constricted_record, columns=boundtime_columns).astype({
        'StartFrame': 'int', 'EndFrame': 'int'})
    df_constrained['time'] = df_constrained['EndFrame'] - df_constrained['StartFrame'] + 1
    df_constrained['type'] = 'C.Diffusion'

    df_strict_bound = pd.DataFrame(bound_strict_record, columns=boundtime_columns).astype({
        'StartFrame': 'int', 'EndFrame': 'int'})
    df_strict_bound['time'] = df_strict_bound['EndFrame'] - df_strict_bound['StartFrame'] + 1
    df_strict_bound['type'] = 'Bound'

    df_fast_diffusion = pd.DataFrame(fastdiff_headers, columns=fastdiff_columns)
    df_fast_diffusion = df_fast_diffusion.dropna(subset=['StartFrame', 'EndFrame'])
    df_fast_diffusion = df_fast_diffusion.astype({'StartFrame': 'int', 'EndFrame': 'int'})
    df_fast_diffusion['time'] = df_fast_diffusion['EndFrame'] - df_fast_diffusion['StartFrame'] + 1
    df_fast_diffusion['type'] = 'FastDiffusion'

    df_all_diffusion = pd.DataFrame(diftime_record, columns=diftime_columns).astype({
        'StartFrame': 'int', 'EndFrame': 'int'})
    df_all_diffusion['time'] = df_all_diffusion['EndFrame'] - df_all_diffusion['StartFrame'] + 1
    df_all_diffusion['type'] = 'SearchTime'

    df_strict_rebind = pd.DataFrame(rebind_strict_time_all, columns=rbtime_columns).astype({
        'StartFrame': 'int', 'EndFrame': 'int'})
    df_strict_rebind['time'] = df_strict_rebind['EndFrame'] - df_strict_rebind['StartFrame'] + 1
    df_strict_rebind['type'] = 'Rebinding'

    combined = pd.concat(
        [df_constrained, df_strict_bound, df_fast_diffusion, df_all_diffusion, df_strict_rebind],
        ignore_index=True,
        sort=False
    )

    # 7) Write the merged table to a single CSV
    out_single_csv = os.path.join(output_path, 'rebind-Events.csv')
    combined.to_csv(out_single_csv, index=False)
    # ───────────────────────────────────────────────────────────────────────
    # ─── Transition Matrix calculation (3×3) ──────────────────────────────────
    transition_matrices = calculate_transition_matrices(tracks)

    # ─── Extract only the 3×3 proportions ───────────────────────────────────────
    prop3x3 = transition_matrices['Proportions']

    # ─── Build the 2×2 “Diffusion vs Bound” proportions ────────────────────────
    abs3x3 = transition_matrices['Absolute Counts']

    from_diff_to_diff   = abs3x3.loc[['From_F.dif','From_C.Dif'], ['To_F.dif','To_C.Dif']].to_numpy().sum()
    from_diff_to_bound  = abs3x3.loc[['From_F.dif','From_C.Dif'], 'To_Bound'].to_numpy().sum()
    from_bound_to_diff  = abs3x3.loc['From_Bound', ['To_F.dif','To_C.Dif']].to_numpy().sum()
    from_bound_to_bound = abs3x3.loc['From_Bound', 'To_Bound']

    counts2x2 = pd.DataFrame(
        [[from_diff_to_diff, from_diff_to_bound],
         [from_bound_to_diff, from_bound_to_bound]],
        index   = ['From Diffusion','From Bound'],
        columns = ['To Diffusion','To Bound']
    )

    row_sums = counts2x2.sum(axis=1)
    props2x2 = counts2x2.div(row_sums, axis=0).fillna(0)

    # ─── Prepare side-by-side printing of 3×3 and 2×2 proportion tables ───────
    str3x3 = prop3x3.to_string()
    str2x2 = props2x2.to_string()

    lines3 = str3x3.split('\n')
    lines2 = str2x2.split('\n')

    width3 = max(len(line) for line in lines3)
    width2 = max(len(line) for line in lines2)
    gap = '   '
    max_lines = max(len(lines3), len(lines2))

    lines3 += [' ' * width3] * (max_lines - len(lines3))
    lines2 += [' ' * width2] * (max_lines - len(lines2))

    combined = [lines3[i] + gap + lines2[i] for i in range(max_lines)]

    # … (code that builds and prints the 3×3/2×2 side-by-side matrix) …

    print_log('\n__Transition Matrix probabilities per Frame__\n')

    # Print “3×3” heading above its three columns, and “2×2” above its two columns
    title3 = '__________________3x3__________________'.center(width3)
    title2 = '__________________2x2__________________'.center(width2)
    print_log(f'{title3}{gap}{title2}')

    for row in combined:
        print_log(row)
    print_log('\n')

    # ─── Now write RESULT_rebind.txt (includes the transition-matrix block above) ───
    with open(log_file) as fin, open(log_result, 'w') as fout:
        active = False
        for line in fin:
            if '[Analysis]' in line:
                active = True
            if active:
                fout.write(line)

    # Also write the CSV‐style results
    with open(csv_result, "w") as csv_file:
        for line in output_result:
            csv_file.write(line)
            csv_file.write('\n')

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

def extract_segments(track, target_label=None, exclude_labels=None):
    segments = []
    current = []
    for spot in track:
        label = int(spot[3])
        if (target_label is not None and label == target_label) or \
           (exclude_labels is not None and label not in exclude_labels):
            current.append(spot)
        elif current:
            segments.append(current)
            current = []
    if current:
        segments.append(current)
    return segments


def extract_segments_from_rebind(track, rebind_durations):
    segments = []
    i = 0
    for duration in rebind_durations:
        while i < len(track) and track[i][3] != 0:
            i += 1
        if i >= len(track):
            break
        segment = []
        while i < len(track) and track[i][3] == 0 and len(segment) < duration:
            segment.append(track[i])
            i += 1
        if len(segment) > 0:
            segments.append(segment)
    return segments


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
        prog='rebind-analysis',
        description='Performs analysis on bound-classified tracks.',
        epilog='Prior scripts: bound_classification.py, gaps_and_fixes.py')
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__ \n\t-> config_path: ' + str(args.config))
    main(args.config)
    print_log("--- %s seconds ---" % (time.time() - start_time))