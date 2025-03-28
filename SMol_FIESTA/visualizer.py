"""
Runnable Script if run as __main__
Generate videos with labelled tracks from bound-decisions or gap-and-fixes.
- Mimic timelapse videos from the tabulated spot positions and their behavior
- Different track behaviors are labelled in different colors.

Input:
    Either of:
        {csv_path}/{output_folder_name}/gaps-and-fixes_decisions.csv
        {csv_path}/{output_folder_name}/bound_decisions.csv
    Determined by parameter: {use_gap_fixed}

Output:
    {csv_path}/{output_folder_name}/#b.tif
        # is replaced with the position of videos in the processing queue

Parameters:
    allowed_track_length_max: max length (frame) of track allowed
    conditional:
        use_gap_fixed: Use tracks processed by gaps_and_fixes.py, if False use tracks only processed by bound_classification.py instead.
"""

import numpy as np
import pandas as pd
from skimage import io as imgio
import os
from natsort import natsorted
import time
import tifffile
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

    outpath = str(os.path.join(csv_path, output_folder_name))

    enable_fixed_particle = False
    use_gap_fixes = configs['toggle']['use_gap_fixed']
    particle_path = 'F:\\MicroscopyTest\\20231210_Dataset\\Fixed_particle\\wt\\particles_result'

    max_frame = configs['track-sorting']['allowed_track_length_max']

    colors = {
        'Cell_Background': [0, 0, 0],
        'Cell_Border': [1, 0, 0],
        'Bound_Center': [102, 204, 51],
        'Bound_Outer': [102, 204, 51],
        'Diffuse_Center': [204, 51, 102],
        'Diffuse_Outer': [204, 51, 102],
        'Constricted_Center': [51, 102, 204],
        'Constricted_Outer': [51, 102, 204],
        'Fixed-Particle': [133, 212, 154],
        'Gap': [0, 0, 0]
    }

    masks = natsorted(get_file_names_with_ext(mask_path, 'png'))
    outlines = natsorted(get_file_names_with_ext(mask_path, 'txt'))
    tracks = pd.read_csv((str(os.path.join(outpath, "bound_decisions.csv"))) if not use_gap_fixes else
                         (str(os.path.join(outpath, "gaps-and-fixes_decisions.csv"))))
    tracks = tracks.loc[:, ~tracks.columns.str.contains('^Unnamed')]

    if enable_fixed_particle:
        particles_files = natsorted(get_file_names_with_ext(particle_path, 'csv'))

    headers = tracks[['Video #', 'Cell', 'Track']].to_numpy()
    tracks = slice_tracks(tracks, headers)
    tracks_by_video = [[] for i in range(len(masks))]
    for track in tracks:
        tracks_by_video[track.iloc[0]['Video #'] - 1].append(track)

    for i in range(len(masks)):
        print('(Video ' + str(i+1) +') -> Mask: ' + masks[i])
        mask = np.swapaxes(imgio.imread(masks[i]), 0, 1)
        videoName= masks[i].split('_')
        if enable_fixed_particle:
            particles = pd.read_csv(particles_files[i]).to_numpy()
        outline = []
        with open(outlines[i], 'r') as file:
            for line in file:
                outline.append(np.array(line.split(',')).astype(int))
        video = inintialize_video(mask, outline, max_frame, colors['Cell_Background'], colors['Cell_Border'], enable_fixed_particle)
        if enable_fixed_particle:
            video = parse_fixed_spots(video, particles, max_frame, i+1, colors['Fixed-Particle'])
        video = parse_tracks(video, tracks_by_video[i], 'Bound', colors)

        video = np.swapaxes(video, 1, 2).astype('uint8')
        save_path = str(os.path.join(outpath, str(i+1) + 'b' + ".tif"))
        print('\t-> Saved to:', save_path)
        tifffile.imwrite(save_path,
                         video, compression='lzw', imagej=True, photometric='rgb', metadata={'axes': 'TYXS', 'mode': 'composite'})
    return

'''
================================================================================================================
KEY FUNCTIONS
================================================================================================================
'''


def parse_tracks(video, tracks_all, key, colors):
    """
    Read all spots and place them in the video with provided position and color based on behavior.
    :param video: (ndarray) A 4d numpy array encoding the color of each pixel over time.
    :param tracks_all: (DataFrame) A track (containing spot information) from other scripts.
    :param key: (str) The column label in tracks_all that indicates spot behavior.
    :param colors: (dict(str:list(int)) The dictionary assigning colors (RGB integers) to each track behavior
    :return:
    """

    for tracks in tracks_all:
        record = (int(tracks.iloc[0]['x']), int(tracks.iloc[0]['y']))
        for iter in range(len(tracks.index)):
            spot = tracks.iloc[iter]
            mark = spot[key]
            frame, x, y = int(spot['Frame']), int(np.round(spot['x'])), int(np.round(spot['y']))
            gap = x == -1 or y == -1
            if(gap):
                x, y = record
            else:
                record = (int(np.round(spot['x'])), int(np.round(spot['y'])))
            if(frame >= video.shape[0]): continue # not really necessary
            outer = spot_index(x, y)
            video[frame][x][y] = (
                colors['Gap'].copy()) if gap else (
                colors['Diffuse_Center'].copy()) if mark == 0 else (
                colors['Constricted_Center'].copy()) if mark == 1 else (
                colors['Bound_Center'].copy())
            for x1, y1 in outer:
                try:
                    video[frame][x1][y1] = (
                        colors['Gap'].copy()) if gap else (
                        colors['Diffuse_Outer'].copy()) if mark == 0 else (
                        colors['Constricted_Outer'].copy()) if mark == 1 else (
                        colors['Bound_Outer'].copy())
                except IndexError:
                    continue
    return video


'''
================================================================================================================
VIDEO UTILITY
================================================================================================================
'''


# Parse known fixed spots in the original timelapse (replisome location, etc.)
# Please ignore, deprecated
def parse_fixed_spots(video, particles, max_frame, vid_index, color_particle):
    for iter in range(len(particles)):
        particle = particles[iter]
        nx, ny, xx, xy = particle[7:11].astype(int).tolist()
        for i in np.arange(nx, xx + 1, 1):
            for j in np.arange(ny, xy + 1, 1):
                video[i][j] = color_particle.copy()
    video = np.repeat(video[np.newaxis, :, :, :], max_frame, axis=0)
    return video

# Set up the video background which is a width by height by length array with rgb color
def inintialize_video(mask, outline, max_frame, cell_color, cell_border_color, use_fixed):
    video = np.empty(shape=(mask.shape[0], mask.shape[1], 3))
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if not mask[i][j] == 0:
                video[i][j] = cell_color.copy()

    for i in range(len(outline)):
        x = outline[i][::2]
        y = outline[i][1::2]
        for k in range(len(x)):
            video[x[k]][y[k]] = cell_border_color
    if not use_fixed: video = np.repeat(video[np.newaxis, :, :, :], max_frame, axis=0)
    return video

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

# Transform spot into a square for better clarity
# spot = loc +- 1
def spot_index(x, y):
    return [
        (x+1, y+1),
        (x+1, y),
        (x+1, y-1),
        (x, y+1),
        (x, y-1),
        (x-1, y+1),
        (x-1, y),
        (x-1, y-1)
    ]

'''
================================================================================================================
I/O
================================================================================================================
'''
def get_file_names_with_ext(path:str, ext:str):
    flist = []
    for root, dirs, files in os.walk(path):
        for file in files:
            fname = file.split('.')
            if(fname[-1] == ext):
                flist.append(root + '\\' +  file)
    return flist

'''
================================================================================================================
START
================================================================================================================
'''

# Start Script
if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(
        prog='visualizer',
        description='Generate videos with labelled tracks from bound-decisions or gap-and-fixes',
        epilog='Prior scripts: bound_classification.py, gaps_and_fixes.py')
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__ \n\t-> config_path: ' + str(args.config))
    main(args.config)
    print("--- %s seconds ---" % (time.time() - start_time))