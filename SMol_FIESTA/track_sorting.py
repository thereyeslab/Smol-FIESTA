"""
Runnable Script if run as __main__
Use the given parameters in "script-config.toml" to parse TrackMate output csv.
- Split track by length and gap.
- Apply preliminary filtering by length, gap and overlapping frames.
- Calculate and record distances to previous and later spots in the same track.
- Make a copy of the configuration file.

Input:
    {csv_path}/*_Cell_*_spotsAll.csv
    {mask_path}/*.png

Output:
    {csv_path}/{output_folder_name}/tracks.csv
    {csv_path}/{output_folder_name}/script-config.toml

Parameters:
    allowed_gap_max: max length (frame) allowed per gap, break track if exceeds.
    allowed_track_length_min: min length (frame) of track allowed
    allowed_track_length_max: max length (frame) of track allowed
    dist_range: how many frames before/after to consider for distance calculations
    allowed_overlap: allowed total repeated frames in tracks in a cell, eliminate all tracks if exceeds
    concurrent_max: max concurrent frames in *different* tracks in a cell, eliminate the frames if exceeds
    conditional:
        bacteria_analysis: the analysis is run on small cells, additional filter conditions apply
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
    mask_path = configs['path']['mask_path']
    one_csv_per_cell = configs['path']['one_csv_per_cell']

    allowed_gap_max = configs['track-sorting'].get('allowed_gap_max', 'None')
    allowed_track_length_min = configs['track-sorting'].get('allowed_track_length_min', 'None')
    allowed_track_length_max = configs['track-sorting'].get('allowed_track_length_max', 'None')
    dist_range = configs['track-sorting'].get('dist_range', 5)
    allowed_overlap = configs['track-sorting'].get('allowed_overlap', 'None')
    concurrent_max = configs['track-sorting'].get('concurrent_max', 'None')
    allowed_concurrent = configs['track-sorting'].get('allowed_concurrent', 'None')
    bacteria_analysis = configs['toggle'].get('bacteria_analysis', True)

    dist_none = float('inf')


    # Output Format
    final_list_track_spots = []
    final_list_track_spots_columns = [
        'Video #', 'Video Name', 'Cell', 'Track', 'Frame', 'x', 'y', 'Intensity'
    ]
    dist_index = list(np.arange(-1*dist_range, dist_range + 1, 1))
    del dist_index[dist_range]
    dist_columns = np.array(dist_index).astype(str)
    for i in range(len(dist_columns)):
        dist_columns[i] = 'Dist ' + dist_columns[i]
    dist_columns = list(dist_columns)
    final_list_track_spots_columns += dist_columns

    output_path = str(os.path.join(csv_path, configs['path']['output_folder_name']))
    try:
        shutil.rmtree(output_path)
        os.mkdir(output_path)
    except:
        os.mkdir(output_path)
    logging_setup(output_path, 'track-sorting')

    # config file save
    shutil.copy(config_path, str(os.path.join(output_path + 'script-config.toml')))
    # getting the masks and csv files full paths
    masks = natsorted(get_file_names_with_ext(mask_path, 'png')) # list of paths to masks
    if one_csv_per_cell:
        csv_sorted = csv_name_sort_suffix(csv_path, 'spots') # Dictionary where: Keys are video names (like "video1"),Values are lists of CSV file paths corresponding to that video (one CSV per cell).
        csv_keys = natsorted(list(csv_sorted.keys())) # list of video names (like "video1")
        if not len(csv_sorted.keys()) == len(masks):
            raise ValueError('Different number of Masks and Videos, you may have a video without tracks')
    else:
        spot_csv_combined = natsorted([
        os.path.join(csv_path, f) for f in os.listdir(csv_path)
        if f.endswith('_spots*.csv') and 'cell' not in f.lower()])
        if not len(spot_csv_combined) == len(masks):
            raise ValueError('Different number of Masks and Videos, you may have a video without tracks')

    for i in range(len(masks)):
        print_log('Processing:', masks[i])
        mask = np.swapaxes(imgio.imread(masks[i]), 0, 1)
        n_cell = np.max(mask)
        if one_csv_per_cell:
            video_name = csv_keys[i]
            spot_csv_paths = natsorted(csv_sorted[video_name])  # List of per-cell spot files
            spots_video = index_format(spot_csv_paths, n_cell)  # list of "CSV file paths" for all cells in the video, sorted by cell index
        else:
            spot_csv_path = spot_csv_combined[i]  # one file per video
            video_name = os.path.splitext(os.path.basename(spot_csv_path))[0].replace('_spots', '')
            spots_video = parse_combined_spots_by_mask(mask, spot_csv_path, n_cell) # list of "spot data", each element is a "list" containing spot information for a cell in the current video (sorted by cell index)

        print_log('\t# Cells in Video:', len(spots_video))
        for j in range(n_cell):
            print_log('\t-> Cell', j, end='')
            try:
                if one_csv_per_cell:
                    spots_cell, _ = parse_csv_by_mask(mask, spots_video[j], j + 1) # list, spot data for the current cell in the current video
                else:
                    spots_cell = spots_video[j]  # list, spot data for the current cell in the current video
                    _ = 0
            except Exception as e:
                print_log(f"Error processing CSV file for cell {j}: {e}")
                continue

            if _ is None:
                print_log(' [ NO LIFETIME TRACK ]')
                continue
            print_log(' [ ELIMINATED BY CELL MASK:', _, ']')

            if len(spots_cell) == 0:
                print_log('\t\t [ ALL SPOTS ELIMINATED ]')
                continue

            tracks_cell = track_separation(spots_cell) # list of tracks, each element is a list of spot information for that track
            print_log('\t\t: # Tracks in Cell:', len(tracks_cell))

            if bacteria_analysis:
                # Eliminate repeated frames in the same track <- allowed_overlap parameter
                _ = 0
                for k in range(len(tracks_cell)):
                    tracks_cell[k], __ = eliminate_repeated_frames(tracks_cell[k])
                    _ += __
                print_log('\t\t: # Repeated Spots eliminated:', _)
                if isinstance(allowed_overlap, str):
                    allowed_overlap = 0

                if _ > allowed_overlap:
                    print_log('\t\t: * Overlap within track exceeds threshold: Aborting Cell')
                    continue

                total_concurrent, _ = tabulate_frame_count(tracks_cell)
                if isinstance(allowed_concurrent, str):
                    allowed_concurrent = 0
                if total_concurrent > allowed_concurrent:
                    print_log('\t\t: * Overlap within track exceeds threshold: Aborting Cell')
                    continue

                # Eliminate repeated frames in different tracks if length exceeds <- concurrent_max parameter
                removal_queue, frame_counts = concurrent_count(tracks_cell, concurrent_max)
                if len(removal_queue) > 0:
                    for k in range(len(tracks_cell)):
                        tracks_cell[k] = remove_frames(tracks_cell[k], removal_queue)
                    print_log('\t\t: # Concurrent Frames removed:', len(removal_queue))

            tracks_ind = []
            _ = 0
            _1 = 0
            for track in tracks_cell:
                result = track_splitting_filtered(track, allowed_gap_max, allowed_track_length_min,
                                                  allowed_track_length_max)
                if result is None:
                    print_log("track_splitting_filtered returned None. Skipping this track.")
                    continue
                track_ind, __, __1 = result
                tracks_ind += track_ind
                _ += __
                _1 += __1
            print_log('\t\t: # Splitting:', _, '# Filtered:', _1)
            print_log('\t\t: # Continuous, Individual Tracks:', len(tracks_ind))

            info = [i + 1, video_name, j + 1]
            for k in range(len(tracks_ind)):
                tracks_ind[k] = track_distance_tabulate(tracks_ind[k], dist_index, dist_none)
                for spot in tracks_ind[k]:
                    entry = info + [k + 1] + spot
                    final_list_track_spots.append(entry)

    # Output
    print_log('Saving to csv:', str(os.path.join(output_path, "Intermediates", 'tracks.csv')))
    final_list_track_spots = pd.DataFrame(final_list_track_spots, columns=final_list_track_spots_columns)
    final_list_track_spots.to_csv(str(os.path.join(output_path, "Intermediates", 'tracks.csv')))
    return

'''
================================================================================================================
KEY FUNCTIONS
================================================================================================================
'''

def track_splitting_filtered(track, gap_max, len_min, len_max):
    """
    Track splitting by the allowed maximum frame gap, filter tracks by allowed maximum and minimum length.
    :param track: (list(list(float)) A track (list of spot information) identified by Trackmate.
    :param gap_max: (int) allowed maximum frame gap before splitting track
    :param len_min: (int) allowed minimum frame length for a track
    :param len_max: (int) allowed maximum frame length for a track
    :return: A list of split tracks, and counts of filtering action for logging.
    """

    res = []
    count_split = 0
    count_filter = 0
    record = []
    if isinstance(len_min, str):
        len_min = 0
    if isinstance(len_max, str):
        len_max = 1000000
    if isinstance(gap_max, str):
        gap_max = 1000000

    # Split
    for i in range(len(track) - 1):
        record.append(track[i])
        if track[i+1][0] - track[i][0] > gap_max:
            res.append(record.copy())
            count_split += 1
            record = []
    try:
        record.append(track[len(track) - 1])
    except IndexError:
        # Handle the IndexError here (e.g., print a message or take corrective action)
        print("IndexError: list index out of range occurred. Skipping this track.")
        return None
    except TypeError:
        # Handle the TypeError here
        print("TypeError: cannot unpack non-iterable NoneType object. Skipping this track.")
        return None

    res.append(record.copy())

    # Filter by track length
    i = 0
    while i < len(res):
        if(len(res[i]) < len_min or len(res[i]) > len_max):
            del res[i]
            count_filter += 1
        else:
            i += 1

    return res, count_split, count_filter

'''
================================================================================================================
TRACKS
================================================================================================================
'''

def distance(p1, p2):
    return np.sqrt(np.power(p1[0] - p2[0], 2) + np.power(p1[1] - p2[1], 2))

# count total number of concurrent track per frame, update frame count dict also
def tabulate_frame_count(tracks):
    total = 0
    frames = {}
    for track in tracks:
        for spot in track:
            if spot[0] in frames:
                if frames[spot[0]] == 1: total += 1
                frames[spot[0]] += 1
            else:
                frames[spot[0]] = 1
    return total, frames


# Remove concurrent spots based on the length of overlap
def concurrent_count(tracks, max_concurrent):
    if isinstance(max_concurrent, str):
        max_concurrent = 1000000
    frame_counts = {}
    for i in range(len(tracks)):
        track = tracks[i]
        for f in range(len(track)):
            frame = track[f][0]
            if not frame in frame_counts:
                frame_counts[frame] = 1
            else:
                frame_counts[frame] += 1

    removal = []
    que = []
    frame_keys = sorted(list(frame_counts.keys()))
    i = 0
    while i < len(frame_keys):
        count = frame_counts[frame_keys[i]]
        if count > 1:
            que.append(frame_keys[i])
        else:
            if len(que) > 0:
                interval = que[-1] - que[0] + 1
                if interval > max_concurrent:
                    removal += que
                que = []
        i += 1
    if len(que) > 0:
        interval = que[-1] - que[0] + 1
        if interval > max_concurrent:
            removal += que
    return removal, frame_counts

# Removing frames in tracks within a queue
def remove_frames(track, queue):
    res = []
    for f in range(len(track)):
        if track[f][0] in queue:
            continue
        else:
            res.append(track[f])
    return res

# Distance calculations appended to the end of spots
def track_distance_tabulate(track, indices, dist_none):
    for i in range(len(track)):
        for j in range(len(indices)):
            if i + indices[j] < 0 or i + indices[j] >= len(track):
                track[i].append(dist_none)
            else:
                track[i].append(distance(
                    track[i][1:3], track[i + indices[j]][1:3]
                ))
    return track

# Eliminate repeat spots in the same frame within a track
def eliminate_repeated_frames(track):
    count = 0
    res = []
    repeats = []
    scan = 0
    while scan < len(track) - 1:
        frame = track[scan][0]
        if not track[scan + 1][0] == frame:
            if len(repeats) > 0:
                repeats.append(track[scan])
                res.append(decide_spots_elimination(repeats, res[-1] if len(res) > 0 else None, track[scan + 1]))
                count += len(repeats) - 1
                repeats = []
            else:
                res.append(track[scan])
            scan += 1
        else:
            repeats.append(track[scan])
            scan += 1
    if len(repeats) > 0:
        repeats.append(track[scan])
        res.append(decide_spots_elimination(repeats, res[-1] if len(res) > 0 else None, None))
        count += len(repeats) - 1
    return res, count


# Decides which spot to eliminate in the same frame within a track based on their distances
def decide_spots_elimination(repeats, prev, nxt):
    best_index = -1
    best_dist = float('inf')
    for i in range(len(repeats)):
        if prev == None:
            dist = distance(repeats[i][1:3], nxt[1:3])
        elif nxt == None:
            dist = distance(repeats[i][1:3], prev[1:3])

        else:
            dist = distance(repeats[i][1:3], prev[1:3])

        if dist < best_dist:
            best_index = i
            best_dist = dist
    return repeats[best_index]

# Separate each track and sort by frame number
def track_separation(spots):
    n_tracks = int(np.max(np.array(spots)[:, 0]) + 1)
    tracks = []
    for i in range(n_tracks): tracks.append([])
    for i in range(len(spots)):
        tracks[spots[i][0]].append(spots[i])
    i = 0
    while i < len(tracks):
        if len(tracks[i]) == 0:
            del tracks[i]
        else:
            i += 1
    for i in range(len(tracks)):
        for j in range(len(tracks[i])):
            tracks[i][j] = tracks[i][j][1:]
    for i in range(len(tracks)):
        tracks[i].sort(key=lambda x: x[0])
    return tracks

'''
================================================================================================================
MASKS
================================================================================================================
'''

# Read CSV and compare each spot to the mask, eliminate if outside specified cell
def parse_csv_by_mask(mask, csv, index):
    remove = False

    if csv is None:
        return None, None
    try:
        # peek at first line to see if it has any letters
        with open(csv, 'r') as _f:
            first = _f.readline()
        skiprows = 1 if any(c.isalpha() for c in first) else 0

        data = np.loadtxt(csv,
                          delimiter=',',
                          dtype=float,
                          skiprows=skiprows)

    except ValueError as e:
        print(f"Error loading CSV file {csv}: {e}")
        return None, None

    if data.ndim == 1:
        data = np.array([data])

    res = []

    for i in range(data.shape[0]):
        try:
            x = int(round(data[i][2]))
            y = int(round(data[i][3]))
        except (IndexError, ValueError) as e:
            print(f"Error processing row {i} in CSV file {csv}: {e}")
            continue
        if remove:
            try:
                cell = mask[x, y]
                if cell == index:
                    res.append([
                        int(data[i][0]), int(data[i][1]), float(data[i][2]), float(data[i][3]), int(data[i][4])
                    ])

            except (IndexError, ValueError) as e:
                print(f"Error matching cell for row {i} in CSV file {csv}: {e}")
                continue
        else:
            res.append([
                int(data[i][0]), int(data[i][1]), float(data[i][2]), float(data[i][3]), int(data[i][4])
            ])
    return res, data.shape[0] - len(res)

# assign each track file the cell number
def index_format(files, max):
    res = [None]*max
    for file in files:
        index = index_find(file)
        try:
            if(not index == -1):
                res[index - 1] = file
        except:
            raise ValueError('Track csv is for a cell index exceeding number of cells in the mask.')
    return res

# read cell number indicated on file name
def index_find(name):
    info = name.split('_')
    try:
        i = info.index('Cell')
    except ValueError:
        raise ValueError('Track csv file name not formatted correctly: missing \"Cell\" in name.')
    if(i + 1 < len(info)):
        return int(info[i+1])
    else:
        raise ValueError('Track csv file name not formatted correctly: no cell index found.')

def parse_combined_spots_by_mask(mask, csv_path, n_cells):
    """
    Parses a single spots CSV containing data for all cells in the video,
    and assigns each spot to the correct cell using the mask.

    Returns a list, list of "spot data", each element is a "list" containing spot information for a cell in the current video: [list of spots for cell 1, list of spots for cell 2, ..., list of spots for cell n].
    Each spot is represented as a list: [track ID, frame, x, y, intensity].
    """
    try:
        # peek at first line to see if it has any letters
        with open(csv_path, 'r') as _f:
            first = _f.readline()
        skiprows = 1 if any(c.isalpha() for c in first) else 0

        data = np.loadtxt(csv_path,
                          delimiter=',',
                          dtype=float,
                          skiprows=skiprows)

    except ValueError as e:
        print(f"Error loading combined spots CSV {csv_path}: {e}")
        return [None] * n_cells

    if data.ndim == 1:
        data = np.array([data])

    cell_spots = [[] for _ in range(n_cells)]

    for i in range(data.shape[0]):
        try:
            x = int(round(data[i][2]))
            y = int(round(data[i][3]))
            cell_id = int(mask[x, y])
            if 1 <= cell_id <= n_cells:
                spot = [
                    int(data[i][0]),   # track ID
                    int(data[i][1]),   # frame
                    float(data[i][2]), # x
                    float(data[i][3]), # y
                    int(data[i][4])    # intensity
                ]
                cell_spots[cell_id - 1].append(spot)
        except (IndexError, ValueError) as e:
            print(f"Error processing spot {i} in {csv_path}: {e}")
            continue

    return cell_spots

'''
================================================================================================================
FILE HANDLING
================================================================================================================
'''

# Filter csv with suffix -> "spotsAll" by default
def csv_name_sort_suffix(path: str, suffix:str='spotsAll') -> dict:
    flist = get_file_names_with_ext(path, 'csv')
    csv_sorted = {}
    for file in flist:
        fname = str(os.path.splitext(os.path.basename(file))[0]).split('_')
        if len(fname) < 4:
            continue
        if 'Cell' not in fname:
            continue
        ind = fname.index('Cell')
        video = str('_').join(fname[:ind])
        if (not video in csv_sorted):
            csv_sorted[video] = []
        if suffix in fname[ind + 2]:
            csv_sorted[video].append(file)
    return csv_sorted

# get all files under path with extension
def get_file_names_with_ext(path: str, ext: str):
    """
    Recursively collect all files under `path` with extension `ext`,
    ignoring any hidden files or directories (names starting with '.').
    """
    flist = []
    for root, dirs, files in os.walk(path):
        # drop hidden dirs
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        for fname in files:
            # skip hidden files
            if fname.startswith('.'):
                continue
            if fname.lower().endswith(f".{ext.lower()}"):
                flist.append(os.path.join(root, fname))
    return flist



'''
================================================================================================================
START
================================================================================================================
'''


# Setup Logging
def logging_setup(path:str, script_name:str):
    logs_dir = os.path.join(path, 'Logs')
    Inter = os.path.join(path, 'Intermediates')
    os.makedirs(Inter, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    log_file = str(os.path.join(logs_dir, 'LOG_' + script_name + '.txt'))
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
        prog='track-sorting',
        description='Use the given parameters in "script-config.toml" to parse TrackMate output csv.',
        epilog='')
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__ \n\t-> config_path: ' + str(args.config))
    main(args.config)
    print_log("--- %s seconds ---" % (time.time() - start_time))