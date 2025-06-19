import sys
import threading
from ij import IJ
from ij import WindowManager
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import ThresholdDetectorFactory
from fiji.plugin.trackmate.detection import DogDetectorFactory
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.util import TMUtils
from ij.plugin.frame import RoiManager
from ij.gui import PolygonRoi
from ij.gui import Roi
from ij.plugin import RoiEnlarger
from java.awt import FileDialog
import csv
import os
import math
import time
import itertools as itto
import shutil
import ConfigParser as configparser  # Top of file already

# Avoid errors with UTF8 chars generated in TrackMate that will mess with our Fiji Jython.
reload(sys)


def main():
    config_file = os.path.abspath(r"C:\Users\JpRas\OneDrive\Escritorio\Analysis_PIPELINE\new\Trackmate_configFile.txt")
    link_parameters, paths, roi_params, settings, filters, config_path = load_config(config_file)

    files = [f for f in os.listdir(paths['dir_tiff']) if f.endswith(".tif")]
    outlines = []

    if paths['segmentation_enabled']:
        outlines = [f for f in os.listdir(paths['dir_outline']) if f.endswith(".txt")]
        if not len(outlines) == len(files):
            print("File mismatch, aborted")
            return
    else:
        outlines = ['' for _ in files]

    save_analysis_dir = os.path.join(paths['dir_tiff'], 'AnalysisSMolFiesta')
    os.mkdir(save_analysis_dir)

    shutil.copy2(config_file, os.path.join(save_analysis_dir, "used_config.txt"))

    link_parameters['files'] = sorted(files)
    link_parameters['save_analysis_dir'] = save_analysis_dir

    for z in range(len(files)):
        imp = IJ.openImage(os.path.join(paths['dir_tiff'], files[z]))
        tracking(imp, link_parameters, z,
                 run_mode=1,
                 radius=settings['radius'],
                 max_frame_gap=settings['max_frame_gap'],
                 number_threads=settings['number_threads'],
                 filter_set=filters,
                 segmentation_enabled=paths['segmentation_enabled'],
                 outline_file=os.path.join(paths['dir_outline'], sorted(outlines)[z]),
                 expansion_size=roi_params['expansion_size'],
                 roi_ids=roi_params['roi_ids'])
                 
def load_config(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)

    # Convert config sections
    link_parameters = {
        'linking_max_distance': config.getfloat('Tracking', 'linking_max_distance'),
        'gap_closing_max_distance': config.getfloat('Tracking', 'gap_closing_max_distance'),
        'splitting_max_distance': config.getfloat('Tracking', 'splitting_max_distance'),
        'merging_max_distance': config.getfloat('Tracking', 'merging_max_distance'),
        'threshold': config.getfloat('Detection', 'threshold'),
        'alternative_linking_cost_factor': config.getfloat('Tracking', 'alternative_linking_cost_factor'),
        'analysis_suffix': config.get('General', 'analysis_suffix'),
        'detector_type': config.get('Detection', 'detector_type'),
    }

    paths = {
        'dir_tiff': config.get('Paths', 'dir_tiff'),
        'dir_outline': config.get('Paths', 'dir_outline'),
        'segmentation_enabled': config.getboolean('Paths', 'segmentation_enabled'),
    }

    roi_params = {
        'expansion_size': config.getfloat('ROI', 'expansion_size'),
        'roi_ids': config.get('ROI', 'roi_ids'),
    }

    settings = {
        'radius': config.getfloat('Settings', 'radius'),
        'max_frame_gap': config.getint('Settings', 'max_frame_gap'),
        'number_threads': config.getint('Settings', 'number_threads'),
    }

    filters = config.items('Filters') if config.has_section('Filters') else []

    return link_parameters, paths, roi_params, settings, filters, config_path


def tracking(imp, link_parameters, z, run_mode=0, radius=1.5, max_frame_gap=2,
             number_threads=4, filter_set=0, segmentation_enabled=False, outline_file='',
             expansion_size=0.0, roi_ids='all'):
    """
    Performs tracking on the given image with specified parameters.

    Args:
        imp (ImagePlus): The image to be processed.
        link_parameters (dict): Dictionary of tracking parameters.
        z (int): Index of the current file.
        run_mode (int): 0 for initial processing; 1 for filtering steps.
        radius (float): Average radius of spots in pixels.
        max_frame_gap (int): Maximum frame gap for gap closing.
        number_threads (int): Number of threads to use.
        filter_set (int): Selected filter set.
        segmentation_enabled (bool): Whether segmentation is performed.
        outline_file (str): Path to the ROI outline file.
        expansion_size (float): ROI expansion size in pixels.
        roi_ids (str or list): Specific ROI IDs for analysis.
    
    Returns:
        list: A list containing tracking results.
    """
    IJ.run("Fresh Start", "")
    print(imp.getTitle(), outline_file, '*********************************************************************************************')

    if segmentation_enabled:
        rm = load_rois(outline_file, imp, expansion_size, roi_ids)
        cell_rois = rm.getRoisAsArray()
        results = []
        for k, roi in enumerate(cell_rois):
            # Determine ROI ID (using 1-based indexing if roi_ids is "all")
            current_roi_id = k + 1 if roi_ids == "all" else (roi_ids[k] if isinstance(roi_ids, list) else roi_ids)
            imp.setRoi(roi)
            model, settings = configure(imp, link_parameters, radius, max_frame_gap, filter_set)
            settings.setRoi(roi)
            print("Video:", str(z + 1), "Cell:", str(current_roi_id), {roi}, "............................................................")
            stat = run_with_timeout(track_process, 60, model, settings, imp, link_parameters, z, run_mode,
                                    filter_set, number_threads, segmentation_enabled, current_roi_id)
            if stat is not None:
                results.append(stat)
        return results
    else:
        model, settings = configure(imp, link_parameters, radius, max_frame_gap, filter_set)
        return track_process(model, settings, imp, link_parameters, z, run_mode,
                             filter_set, number_threads)


def run_with_timeout(func, timeout, *args, **kwargs):
    """
    Runs a function with a timeout. If it does not complete within the timeout,
    execution moves to the next step.

    Args:
        func (function): The function to run.
        timeout (int): Timeout in seconds.

    Returns:
        The result of the function if completed in time; otherwise, None.
    """
    class InterruptableThread(threading.Thread):
        def __init__(self):
            threading.Thread.__init__(self)
            self.result = None

        def run(self):
            self.result = func(*args, **kwargs)

    thread_instance = InterruptableThread()
    thread_instance.start()
    thread_instance.join(timeout)
    if thread_instance.is_alive():
        print("Timeout: Function {} exceeded {} seconds.".format(func.__name__, timeout))
        return None
    return thread_instance.result


def configure(imp, link_parameters, radius, max_frame_gap, filter_set):
    """
    Configures settings for the tracking process.

    Args:
        imp (ImagePlus): The image to process.
        link_parameters (dict): Dictionary of tracking parameters.
        radius (float): Average radius of spots.
        max_frame_gap (int): Maximum frame gap for gap closing.
        filter_set (int): Selected filter set.
    
    Returns:
        tuple: (model, settings)
    """
    # Create model and set logger
    model = Model()
    model.setLogger(Logger.IJ_LOGGER)

    # Prepare settings object
    settings = Settings(imp)

    # Configure detector based on detector_type in link_parameters
    detector_type = link_parameters['detector_type']
    if detector_type == 'Thresh':
        settings.detectorFactory = ThresholdDetectorFactory()
        settings.detectorSettings = {
            'TARGET_CHANNEL': 1,
            'INTENSITY_THRESHOLD': link_parameters['threshold'],
            'SIMPLIFY_CONTOURS': True,
        }
    elif detector_type == 'Dog':
        settings.detectorFactory = DogDetectorFactory()
        settings.detectorSettings = {
            'DO_SUBPIXEL_LOCALIZATION': True,
            'RADIUS': radius,
            'TARGET_CHANNEL': 1,
            'THRESHOLD': link_parameters['threshold'],
            'DO_MEDIAN_FILTERING': False,
        }
    elif detector_type == 'Log':
        settings.detectorFactory = LogDetectorFactory()
        settings.detectorSettings = {
            'DO_SUBPIXEL_LOCALIZATION': True,
            'RADIUS': radius,
            'TARGET_CHANNEL': 1,
            'THRESHOLD': link_parameters['threshold'],
            'DO_MEDIAN_FILTERING': False,
        }

    # Configure tracker settings
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = link_parameters['linking_max_distance']
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = link_parameters['gap_closing_max_distance']
    settings.trackerSettings['SPLITTING_MAX_DISTANCE'] = link_parameters['splitting_max_distance']
    settings.trackerSettings['MERGING_MAX_DISTANCE'] = link_parameters['merging_max_distance']
    settings.trackerSettings['MAX_FRAME_GAP'] = max_frame_gap
    settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
    settings.trackerSettings['ALLOW_TRACK_MERGING'] = False
    settings.trackerSettings['ALTERNATIVE_LINKING_COST_FACTOR'] = link_parameters['alternative_linking_cost_factor']

    # Add analyzers and an initial spot filter
    settings.addAllAnalyzers()
    settings.initialSpotFilterValue = 1.0

    def parse_filter_string(filter_str):
        parts = [x.strip() for x in filter_str.split(',')]
        return FeatureFilter(parts[0], float(parts[1]), parts[2].lower() == 'true')

    for name, val in filter_set:
        if name.startswith('spot_filter'):
            settings.addSpotFilter(parse_filter_string(val))
        elif name.startswith('track_filter'):
            settings.addTrackFilter(parse_filter_string(val))

    return model, settings



def track_process(model, settings, imp, link_parameters, z, run_mode, filter_set,
                  number_threads, segmentation_enabled=False, roi_id=0):
    """
    Processes tracking with the given model and settings.

    Args:
        model (Model): The TrackMate model.
        settings (Settings): The TrackMate settings.
        imp (ImagePlus): The image to be processed.
        link_parameters (dict): Dictionary of tracking parameters.
        z (int): Index of the current file.
        run_mode (int): 0 for initial processing; 1 for filtering.
        filter_set (int): Selected filter set.
        number_threads (int): Number of threads to use.
        segmentation_enabled (bool): Whether segmentation is enabled.
        roi_id (int): The ID of the current ROI.
    
    Returns:
        None
    """
    trackmate = TrackMate(model, settings)
    trackmate.setNumThreads(number_threads)

    # Execute checks and processing
    if not trackmate.checkInput():
        print('Input error in TrackMate configuration')
        return
    if not trackmate.process():
        print('Processing error in TrackMate')
        return

    spt_m = model.getSpots()
    tracks_found = model.getTrackModel().trackIDs(True)
    tot_spts = spt_m.getNSpots(True)
    print("Spot total:", tot_spts, "Tracks total:", len(tracks_found))
    model.getLogger().log('Found ' + str(model.getTrackModel().nTracks(True)) + ' tracks.')

    sm = SelectionModel(model)
    ds = DisplaySettingsIO.readUserDefault()

    # Uncomment the following lines to display the viewer if needed.
    displayer = HyperStackDisplayer(model, sm, imp, ds)
    displayer.render()

    fm = model.getFeatureModel()

    # Initialize lists for storing track features
    Track_IDs = []
    mean_sp = []
    med_sp = []
    min_sp = []
    max_sp = []
    std_sp = []
    mean_q = []
    med_q_tr = []
    min_q_tr = []
    max_q_tr = []
    std_q_tr = []
    x_lc = []
    y_lc = []
    mn_int = []
    inten = []
    tr_dur = []
    tr_start = []
    tr_fin = []
    spt_tr = []
    spt_widt = []
    x_tr = []
    y_tr = []
    tr_identifi = []
    tr_fram = []
    poolIDs = []

    # Iterate over visible tracks
    for track_id in model.getTrackModel().trackIDs(True):
        v = fm.getTrackFeature(track_id, 'TRACK_MEAN_SPEED')
        med_v = fm.getTrackFeature(track_id, 'TRACK_MEDIAN_SPEED')
        min_v = fm.getTrackFeature(track_id, 'TRACK_MIN_SPEED')
        max_v = fm.getTrackFeature(track_id, 'TRACK_MAX_SPEED')
        std_v = fm.getTrackFeature(track_id, 'TRACK_STD_SPEED')
        q = fm.getTrackFeature(track_id, 'TRACK_MEAN_QUALITY')
        med_q = fm.getTrackFeature(track_id, 'TRACK_MEDIAN_QUALITY')
        min_q = fm.getTrackFeature(track_id, 'TRACK_MIN_QUALITY')
        max_q = fm.getTrackFeature(track_id, 'TRACK_MAX_QUALITY')
        std_q = fm.getTrackFeature(track_id, 'TRACK_STD_QUALITY')
        dura = fm.getTrackFeature(track_id, 'TRACK_DURATION')
        start_tr = fm.getTrackFeature(track_id, 'TRACK_START')
        fin_tr = fm.getTrackFeature(track_id, 'TRACK_STOP')
        spts = fm.getTrackFeature(track_id, 'NUMBER_SPOTS')
        tr_identif = fm.getTrackFeature(track_id, 'TRACK_ID')
        identi = [int(track_id)] * int(spts)
        x_loc = fm.getTrackFeature(track_id, 'TRACK_X_LOCATION')
        y_loc = fm.getTrackFeature(track_id, 'TRACK_Y_LOCATION')
        poolIDs.append(track_id)

        track = model.getTrackModel().trackSpots(track_id)
        track_int = []
        track_x = []
        track_y = []
        fram = []
        qlist = []
        track_shp = []
        for spot in track:
            x = spot.getFeature('POSITION_X')
            y = spot.getFeature('POSITION_Y')
            t = spot.getFeature('FRAME')
            q_val = spot.getFeature('QUALITY')
            tot_int_spot = spot.getFeature('TOTAL_INTENSITY_CH1')
            wid = spot.getFeature('RADIUS')
            track_int.append(tot_int_spot)
            track_x.append(x)
            track_y.append(y)
            fram.append(t)
            qlist.append(q_val)
            track_shp.append(wid)

        Track_IDs.append(track_id)
        mean_sp.append(v)
        med_q_value = qlist[int(round(len(qlist) / 2))]
        min_q_value = min(qlist)
        max_q_value = max(qlist)
        std_q_value = math.sqrt(sum((x - sum(qlist) / len(qlist)) ** 2 for x in qlist) / len(qlist))
        med_sp.append(med_v)
        min_sp.append(min_v)
        max_sp.append(max_v)
        std_sp.append(std_v)
        mean_q.append(q)
        med_q_tr.append(med_q_value)
        min_q_tr.append(min_q_value)
        max_q_tr.append(max_q_value)
        std_q_tr.append(std_q_value)
        x_loc_avg = sum(track_x) / len(track_x)
        y_loc_avg = sum(track_y) / len(track_y)
        x_lc.append(x_loc_avg)
        y_lc.append(y_loc_avg)
        mean_track_int = sum(track_int) / len(track_int)
        mn_int.append(mean_track_int)
        inten.extend(track_int)
        tr_dur.append(dura)
        tr_start.append(start_tr)
        tr_fin.append(fin_tr)
        spt_tr.append(spts)
        spt_widt.append(sum(track_shp) / len(track_shp))
        tr_fram.extend(fram)
        x_tr.extend(track_x)
        y_tr.extend(track_y)
        tr_identifi.extend(identi)

    qTracks = TMUtils.otsuThreshold(mean_q)
    qSpots = TMUtils.otsuThreshold(x_tr)

    # Organize spot-level data
    spt_all_fin = [tr_identifi, tr_fram, x_tr, y_tr, inten]
    spt_all_fin_2 = list(zip(*spt_all_fin))
    data_final = [Track_IDs, spt_tr, spt_widt, mean_sp, max_sp, min_sp, med_sp, std_sp,
                  mean_q, max_q_tr, min_q_tr, med_q_tr, std_q_tr, tr_dur, tr_start, tr_fin, x_lc, y_lc]
    data_final_2 = list(zip(*data_final))

    Count = {item: spt_all_fin[1].count(item) for item in spt_all_fin[1]}
    CountList = [(k, v) for k, v in Count.items() if v > 0.9]
    print('Number of spots:', len(Count), "Concurrent Localizations:", len(CountList))

    if len(model.getTrackModel().trackIDs(True)) > 0 and spt_m.getNSpots(True) < 2000000000000:
        base_filename = link_parameters['files'][z][:-4]
        analysis_suffix = link_parameters['analysis_suffix']
        save_dir = link_parameters['save_analysis_dir']
        tracks_filename = os.path.join(save_dir, "{}_Cell_{}_tracks{}.csv".format(base_filename, roi_id, analysis_suffix))
        spots_filename = os.path.join(save_dir, "{}_Cell_{}_spots{}.csv".format(base_filename, roi_id, analysis_suffix))

        with open(tracks_filename, "w") as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerows(data_final_2)
        with open(spots_filename, "w") as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerows(spt_all_fin_2)


        IJ.log('Success')
        return


def load_rois(file_name, imp, expansion_size, roi_ids):
    """
    Loads ROIs from a text file and expands them.

    Args:
        file_name (str): Path to the outline file.
        imp (ImagePlus): The image to set the ROIs on.
        expansion_size (float): Size for ROI expansion in pixels.
        roi_ids (str or list): Specific ROI IDs for analysis.
    
    Returns:
        RoiManager: ROI manager containing the loaded and expanded ROIs.
    """
    rm = RoiManager.getInstance()
    if not rm:
        rm = RoiManager()
    rm.reset()

    with open(file_name, "r") as textfile:
        xy = [list(map(int, line.rstrip().split(","))) for line in textfile]

    for i in range(len(xy)):
        if roi_ids == "all" or i + 1 in roi_ids:
            x = xy[i][::2]
            y = xy[i][1::2]
            imp.setRoi(PolygonRoi(x, y, Roi.POLYGON))
            roi = imp.getRoi()
            roi = RoiEnlarger.enlarge(roi, expansion_size)
            rm.addRoi(roi)
    return rm



if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- {} seconds ---".format(time.time() - start_time))
