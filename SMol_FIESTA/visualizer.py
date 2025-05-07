#!/usr/bin/env python3
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
        new_arr[::2] -= min_row
        new_arr[1::2] -= min_col
        new_outline.append(new_arr)
    return cropped_mask, (min_row, min_col), new_outline

def draw_spot(image, center, diameter, color):
    """Draws a filled circle (spot) on the image with subpixel localization."""
    radius = diameter / 2.0
    r0 = int(np.floor(center[0] - radius))
    r1 = int(np.ceil(center[0] + radius))
    c0 = int(np.floor(center[1] - radius))
    c1 = int(np.ceil(center[1] + radius))
    H, W, _ = image.shape
    r0, c0 = max(r0, 0), max(c0, 0)
    r1, c1 = min(r1, H), min(c1, W)
    for i in range(r0, r1):
        for j in range(c0, c1):
            d = np.hypot(i + 0.5 - center[0], j + 0.5 - center[1])
            if d <= radius:
                image[i, j, :] = color
    return image

def draw_tail(image, start, end, tail_diam, color):
    """Draws a motion-blur tail from start to end."""
    distance = np.linalg.norm(np.array(end) - np.array(start))
    num_steps = max(int(np.ceil(distance)), 2)
    for s in np.linspace(0, 1, num_steps):
        interp = (1 - s) * np.array(start) + s * np.array(end)
        image = draw_spot(image, interp, tail_diam, color)
    return image

def inintialize_video(mask, outline, max_frame, cell_color, cell_border_color, use_fixed):
    # Do NOT crop: keep full original frame so spots align correctly
    cropped_mask = mask
    crop_offset = (0, 0)
    new_outline = outline

    H, W = cropped_mask.shape
    video_bg = np.zeros((H, W, 3), dtype=np.uint8)
    for i in range(H):
        for j in range(W):
            if cropped_mask[i, j] != 0:
                video_bg[i, j] = cell_color
    for arr in new_outline:
        xs = arr[::2].astype(int)
        ys = arr[1::2].astype(int)
        for x, y in zip(xs, ys):
            if 0 <= x < H and 0 <= y < W:
                video_bg[x, y] = cell_border_color

    if not use_fixed:
        video = np.repeat(video_bg[np.newaxis, :, :, :], max_frame, axis=0)
    else:
        video = video_bg[np.newaxis, :, :, :]
    return video, crop_offset

def crop_video_to_spot_bbox(video, tracks_for_video, crop_offset, margin=0):
    # No-op now that we keep full frame; return video unchanged
    return video

# ------------------------- Updated parse_tracks -------------------------

def parse_tracks(video, tracks_all, key, colors, crop_offset, spot_diameter):
    tail_effect_enabled = configs.get('visualizer', {}).get('tail_effect_enabled', True)
    tail_diam = spot_diameter / 4

    for tracks in tracks_all:
        prev_adj = None
        for _, spot in tracks.iterrows():
            frame = int(spot['Frame'])
            x, y = spot['x'], spot['y']
            mark = spot[key]
            gap = (x == -1 or y == -1)

            if gap:
                if prev_adj is None:
                    continue
                adj_x, adj_y = prev_adj
            else:
                adj_x = x - crop_offset[0]
                adj_y = y - crop_offset[1]

            if frame >= video.shape[0]:
                continue

            if gap:
                color_used = colors['Gap']
            else:
                if mark == 0:
                    color_used = colors['Diffuse']
                elif mark == 1:
                    color_used = colors['Constricted']
                else:
                    color_used = colors['Bound']

            if (not gap) and prev_adj is not None and mark in (0, 1) and tail_effect_enabled:
                tail_end = (adj_x, adj_y)
                video[frame] = draw_tail(video[frame], prev_adj, tail_end, tail_diam, color_used)

            video[frame] = draw_spot(video[frame], (adj_x, adj_y), spot_diameter, color_used)
            if not gap:
                prev_adj = (adj_x, adj_y)

    return video

# ------------------------- Slice and File Helpers -------------------------

def slice_tracks(tracks, headers):
    indices = []
    last = np.array([-1, -1, -1])
    for i, row in enumerate(headers):
        if not np.array_equal(row, last):
            indices.append(i)
            last = row.copy()
    indices.append(len(headers))
    return [tracks.iloc[indices[i]:indices[i+1]] for i in range(len(indices)-1)]

def get_file_names_with_ext(path: str, ext: str):
    flist = []
    for root, _, files in os.walk(path):
        for file in files:
            if file.lower().endswith(f'.{ext.lower()}'):
                flist.append(os.path.join(root, file))
    return flist

# ------------------------- Main -------------------------

def main(config_path: str = None):
    global configs
    if not config_path:
        base = os.path.dirname(__file__)
        config_path = os.path.join(base, 'script-config.toml')
    with open(config_path, 'rb') as f:
        configs = tomllib.load(f)

    csv_path = configs['path']['csv_path']
    mask_path = configs['path']['mask_path']
    output_folder_name = configs['path']['output_folder_name']
    data_dir = os.path.join(csv_path, output_folder_name)
    outpath = os.path.join(data_dir, 'Visualizer')
    os.makedirs(outpath, exist_ok=True)

    use_gap_fixes = configs['toggle']['use_gap_fixed']
    max_frame = configs['track-sorting']['allowed_track_length_max']

    colors = {
        'Cell_Background': [0, 0, 0],
        'Cell_Border': [86, 86, 86],
        'Bound': [102, 204,  51],
        'Diffuse': [204,  51, 102],
        'Constricted': [ 51, 102, 204],
        'Gap': [  0,   0,   0],
    }


    masks = natsorted(get_file_names_with_ext(mask_path, 'png'))
    outlines = natsorted(get_file_names_with_ext(mask_path, 'txt'))
    input_file = 'gaps-and-fixes_decisions.csv' if use_gap_fixes else 'bound_decisions.csv'
    tracks_df = pd.read_csv(os.path.join(data_dir, input_file))

    headers = tracks_df[['Video #', 'Cell', 'Track']].to_numpy()
    sliced = slice_tracks(tracks_df, headers)
    tracks_by_video = [[] for _ in masks]
    for tr in sliced:
        vid = int(tr.iloc[0]['Video #']) - 1
        tracks_by_video[vid].append(tr)

    spot_diameter = configs.get('visualizer', {}).get('spot_diameter', 6.0)

    for i, mask_file in enumerate(masks):
        print(f'(Video {i+1}) -> Mask: {mask_file}')
        mask = np.swapaxes(imgio.imread(mask_file), 0, 1)

        outline = []
        with open(outlines[i], 'r') as f:
            for line in f:
                outline.append(np.fromstring(line, sep=',', dtype=int))

        video, crop_offset = inintialize_video(
            mask, outline, max_frame,
            colors['Cell_Background'], colors['Cell_Border'],
            use_fixed=False
        )

        video = parse_tracks(video, tracks_by_video[i], 'Bound', colors, crop_offset, spot_diameter)

        # Truncate to last spot + 5
        if tracks_by_video[i]:
            last_frame = max(tr['Frame'].max() for tr in tracks_by_video[i])
            new_len = min(last_frame + 5, video.shape[0])
            video = video[:new_len]

        # No spatial cropping—keep full frame
        # video = crop_video_to_spot_bbox(...)

        video = np.swapaxes(video, 1, 2).astype('uint8')
        out_name = os.path.splitext(os.path.basename(mask_file))[0] + '.tif'
        save_path = os.path.join(outpath, out_name)

        print(f'\t-> Saved to: {save_path}')
        tifffile.imwrite(
            save_path,
            video,
            compression='lzw',
            imagej=True,
            photometric='rgb',
            metadata={'axes': 'TYXS', 'mode': 'composite'}
        )

if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser(
        prog='visualizer',
        description='Generate videos with labelled tracks from bound-decisions or gap-and-fixes',
        epilog='Prior scripts: bound_classification.py, gaps_and_fixes.py'
    )
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__')
    print(f'\t-> config_path: {args.config}')
    main(args.config)
    print(f"--- {time.time() - start:.2f} seconds ---")
