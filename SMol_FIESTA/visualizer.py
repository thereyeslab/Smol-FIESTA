"""
Runnable Script if run as __main__
Generate videos with labelled tracks from bound-decisions or gap-and-fixes.
- Mimic timelapse videos from the tabulated spot positions and their behavior.
- Spots for each behavior are drawn using a single color per class.
- A motion-blur tail is drawn for moving molecules (classes 0 and 1) to indicate motion.
- Video is truncated to 5 frames after the last spot.
- Output file is named using the mask name.
- Configurable parameter: spot_diameter.

Input:
    Either of:
        {csv_path}/{output_folder_name}/gaps-and-fixes_decisions.csv
        {csv_path}/{output_folder_name}/bound_decisions.csv
    Determined by parameter: {use_gap_fixed}

Output:
    {csv_path}/{output_folder_name}/{mask_name}.tif

Parameters:
    allowed_track_length_max: max length (frame) of track allowed
    conditional:
        use_gap_fixed: Use tracks processed by gaps_and_fixes.py; if False, use tracks only processed by bound_classification.py.
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

# ------------------------- Helper Functions -------------------------
def crop_mask_and_outline(mask, outline):
    nonzero = np.argwhere(mask != 0)
    if nonzero.size == 0:
        return mask, (0, 0), outline
    min_row, min_col = nonzero.min(axis=0)
    max_row, max_col = nonzero.max(axis=0)
    cropped_mask = mask[min_row:max_row+1, min_col:max_col+1]
    new_outline = []
    for arr in outline:
        new_arr = arr.copy()
        new_arr[::2] = new_arr[::2] - min_row
        new_arr[1::2] = new_arr[1::2] - min_col
        new_outline.append(new_arr)
    return cropped_mask, (min_row, min_col), new_outline

def draw_spot(image, center, diameter, color):
    """Draws a filled circle (spot) on the image with subpixel localization using the given diameter."""
    radius = diameter / 2.0
    # Set bounding box (no extra padding)
    r0 = int(np.floor(center[0] - radius))
    r1 = int(np.ceil(center[0] + radius))
    c0 = int(np.floor(center[1] - radius))
    c1 = int(np.ceil(center[1] + radius))
    H, W, _ = image.shape
    r0, c0 = max(r0, 0), max(c0, 0)
    r1, c1 = min(r1, H), min(c1, W)
    for i in range(r0, r1):
        for j in range(c0, c1):
            # Use center offset of 0.5 for subpixel accuracy
            d = np.sqrt((i + 0.5 - center[0])**2 + (j + 0.5 - center[1])**2)
            if d <= radius:
                image[i, j, :] = color
    return image

def draw_tail(image, start, end, tail_diam, color):
    """
    Draws a tail (motion blur effect) from start to end using the same color.
    The tail is drawn as a series of spots along the interpolated path.
    """
    distance = np.linalg.norm(np.array(end) - np.array(start))
    num_steps = max(int(np.ceil(distance)), 2)
    for s in np.linspace(0, 1, num_steps):
        interp = (1 - s) * np.array(start) + s * np.array(end)
        image = draw_spot(image, interp, tail_diam, color)
    return image

def inintialize_video(mask, outline, max_frame, cell_color, cell_border_color, use_fixed):
    # Crop mask and outline to minimal bounding box.
    cropped_mask, crop_offset, new_outline = crop_mask_and_outline(mask, outline)
    H, W = cropped_mask.shape
    video_bg = np.empty((H, W, 3), dtype=np.uint8)
    for i in range(H):
        for j in range(W):
            video_bg[i, j] = cell_color.copy() if cropped_mask[i, j] != 0 else [0, 0, 0]
    for arr in new_outline:
        x = arr[::2]
        y = arr[1::2]
        for k in range(len(x)):
            if 0 <= x[k] < H and 0 <= y[k] < W:
                video_bg[x[k], y[k]] = cell_border_color
    if not use_fixed:
        video = np.repeat(video_bg[np.newaxis, :, :, :], max_frame, axis=0)
    else:
        video = video_bg
    return video, crop_offset

def crop_video_to_spot_bbox(video, tracks_for_video, crop_offset, margin=0):
    all_x, all_y = [], []
    for track in tracks_for_video:
        xs = track['x'].values
        ys = track['y'].values
        valid = (xs != -1) & (ys != -1)
        all_x.extend(xs[valid].tolist())
        all_y.extend(ys[valid].tolist())
    if len(all_x) == 0 or len(all_y) == 0:
        return video
    min_x = max(int(np.floor(min(all_x))) - margin, 0)
    max_x = int(np.ceil(max(all_x))) + margin
    min_y = max(int(np.floor(min(all_y))) - margin, 0)
    max_y = int(np.ceil(max(all_y))) + margin
    crop_min_row = int(min_x - crop_offset[0])
    crop_max_row = int(max_x - crop_offset[0])
    crop_min_col = int(min_y - crop_offset[1])
    crop_max_col = int(max_y - crop_offset[1])
    T, H, W, C = video.shape
    crop_min_row = max(crop_min_row, 0)
    crop_min_col = max(crop_min_col, 0)
    crop_max_row = min(crop_max_row, H)
    crop_max_col = min(crop_max_col, W)
    return video[:, crop_min_row:crop_max_row, crop_min_col:crop_max_col, :]

# ------------------------- Updated parse_tracks -------------------------
def parse_tracks(video, tracks_all, key, colors, crop_offset, spot_diameter):
    """
    Draws spots on the video based on track positions.
    - Uses a single color per behavior.
    - Draws a tail (motion blur) for moving molecules of classes 0 and 1.
    """
    # Tail effect: tail diameter is fixed to 1/4 of spot_diameter.
    tail_effect_enabled = configs.get('visualizer', {}).get('tail_effect_enabled', True)
    tail_diam = spot_diameter / 4

    for tracks in tracks_all:
        prev_adj = None  # previous valid adjusted position
        for idx in range(len(tracks.index)):
            spot = tracks.iloc[idx]
            mark = spot[key]
            frame = int(spot['Frame'])
            x = float(spot['x'])
            y = float(spot['y'])
            gap = (x == -1 or y == -1)
            if gap:
                if prev_adj is None:
                    continue
                adj_x, adj_y = prev_adj
            else:
                # Adjust based on crop_offset (no scaling)
                adj_x = x - crop_offset[0]
                adj_y = y - crop_offset[1]
            if frame >= video.shape[0]:
                continue
            # Select single color per behavior
            if gap:
                color_used = colors['Gap'].copy()
            else:
                if mark == 0:
                    color_used = colors['Diffuse'].copy()
                elif mark == 1:
                    color_used = colors['Constricted'].copy()
                else:
                    color_used = colors['Bound'].copy()
            # Draw tail for classes 0 and 1 if enabled and previous position exists
            if (not gap) and prev_adj is not None and (mark == 0 or mark == 1) and tail_effect_enabled:
                tail_end = (prev_adj[0] + (adj_x - prev_adj[0]), prev_adj[1] + (adj_y - prev_adj[1]))
                video[frame] = draw_tail(video[frame], prev_adj, tail_end, tail_diam, color=color_used)
            # Draw the spot using the single color
            video[frame] = draw_spot(video[frame], (adj_x, adj_y), spot_diameter, color=color_used)
            if not gap:
                prev_adj = (adj_x, adj_y)
    return video

# ------------------------- Slice and File Helpers -------------------------
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
        tracks_sliced.append(tracks.iloc[indices[i]:indices[i+1], :])
    return tracks_sliced

def get_file_names_with_ext(path: str, ext: str):
    flist = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.split('.')[-1] == ext:
                flist.append(os.path.join(root, file))
    return flist

# ------------------------- Main -------------------------
def main(config_path: str = None):
    global configs  # to allow access in parse_tracks
    if not config_path:
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, 'script-config.toml')
    with open(config_path, 'rb') as config_file:
        configs = tomllib.load(config_file)

    # Paths and basic parameters
    csv_path = configs['path']['csv_path']
    mask_path = configs['path']['mask_path']
    output_folder_name = configs['path']['output_folder_name']
    outpath = os.path.join(csv_path, output_folder_name)

    enable_fixed_particle = False
    use_gap_fixes = configs['toggle']['use_gap_fixed']
    particle_path = 'F:\\MicroscopyTest\\20231210_Dataset\\Fixed_particle\\wt\\particles_result'
    max_frame = configs['track-sorting']['allowed_track_length_max']

    # Updated color dictionary (single color per behavior)
    colors = {
        'Cell_Background': [0, 0, 0],
        'Cell_Border': [1, 0, 0],
        'Bound': [102, 204, 51],
        'Diffuse': [204, 51, 102],
        'Constricted': [51, 102, 204],
        'Fixed-Particle': [133, 212, 154],
        'Gap': [0, 0, 0]
    }

    masks = natsorted(get_file_names_with_ext(mask_path, 'png'))
    outlines = natsorted(get_file_names_with_ext(mask_path, 'txt'))
    input_file = "bound_decisions.csv" if not use_gap_fixes else "gaps-and-fixes_decisions.csv"
    tracks = pd.read_csv(os.path.join(outpath, input_file))
    tracks = tracks.loc[:, ~tracks.columns.str.contains('^Unnamed')]

    if enable_fixed_particle:
        particles_files = natsorted(get_file_names_with_ext(particle_path, 'csv'))

    headers = tracks[['Video #', 'Cell', 'Track']].to_numpy()
    tracks = slice_tracks(tracks, headers)
    tracks_by_video = [[] for _ in range(len(masks))]
    for track in tracks:
        tracks_by_video[track.iloc[0]['Video #'] - 1].append(track)

    # Visualization parameter: single spot_diameter (default 3.0)
    spot_diameter = configs.get('visualizer', {}).get('spot_diameter', 6.0)

    for i in range(len(masks)):
        mask_filename = os.path.basename(masks[i])
        print('(Video ' + str(i+1) +') -> Mask: ' + masks[i])
        mask = np.swapaxes(imgio.imread(masks[i]), 0, 1)
        if enable_fixed_particle:
            particles = pd.read_csv(particles_files[i]).to_numpy()
        outline = []
        with open(outlines[i], 'r') as file:
            for line in file:
                outline.append(np.array(line.split(',')).astype(int))
        video, crop_offset = inintialize_video(mask, outline, max_frame, colors['Cell_Background'], colors['Cell_Border'], enable_fixed_particle)
        if enable_fixed_particle:
            # Assuming parse_fixed_spots exists if needed
            video = parse_fixed_spots(video, particles, max_frame, i+1, colors['Fixed-Particle'])
        video = parse_tracks(video, tracks_by_video[i], 'Bound', colors, crop_offset, spot_diameter)

        # Truncate video: if no more track positions, cut video at (last spot frame + 5)
        if len(tracks_by_video[i]) > 0:
            max_track_frame = max(track['Frame'].max() for track in tracks_by_video[i])
            new_frame_count = min(max_track_frame + 5, video.shape[0])
            video = video[:new_frame_count]

        # Crop video spatially to region of interest based on track positions (no margin)
        video = crop_video_to_spot_bbox(video, tracks_by_video[i], crop_offset, margin=0)
        video = np.swapaxes(video, 1, 2).astype('uint8')

        # Rename output file to use the mask name (without extension)
        out_filename = os.path.splitext(mask_filename)[0] + ".tif"
        save_path = os.path.join(outpath, out_filename)
        print('\t-> Saved to:', save_path)
        tifffile.imwrite(save_path,
                          video, compression='lzw', imagej=True, photometric='rgb',
                          metadata={'axes': 'TYXS', 'mode': 'composite'})
    return

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
