"""
Runnable Script if run as __main__
Performs MSD calculations with gaussian mixture model.
- Calculate MSD based on spot distances between frames.
- Fit MSD calculations in a gaussian mixture model.
- Generate plots for presenting each population [bound, constrained diffusion, free diffusion] and their MSD distributions.

Input:
    Either of:
        {csv_path}/{output_folder_name}/gaps-and-fixes_decisions.csv
        {csv_path}/{output_folder_name}/bound_decisions.csv
    Determined by parameter: {use_gap_fixed}

Output:
    {csv_path}/{output_folder_name}/Diffusion_Coefficient_Calculation.csv
    {csv_path}/{output_folder_name}/Diffusion_Coefficient_Plots.pdf

Parameters:
    conditional:
        use_gap_fixed: Use tracks processed by gaps_and_fixes.py, if False use tracks only processed by bound_classification.py instead.
"""

import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from scipy.stats import norm
import tomllib
import argparse

def main(config_path:str = None):
    if not config_path:
        __location__ = os.path.realpath(
            os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, 'script-config.toml')
    with open(config_path, 'rb') as config_file:
        configs = tomllib.load(config_file)

    # Load the CSV file
    csv_path = configs['path']['csv_path']
    output_dir = os.path.join(csv_path, configs['path']['output_folder_name'])

    # Conversion factor from pixels to microns
    pixel_to_micron = 0.130

    # Parameters
    frames_per_subdivision = 10  # Number of frames per subdivision
    time_interval = 0.1  # Time interval in seconds
    bin_scaling_factor = 0.75  # Scaling factor for the number of bins

    # Min and max diffusion for each behaviour

    cutoff_bound = 0.02  # Diffusion coefficient cutoff for bound
    cutoff_cdiffusion = 2.5  # Diffusion coefficient cutoff for constrained diffusion
    cutoff_fdiffusion = 5  # Diffusion coefficient cutoff for fast diffusion
    min_bound = 0.0001  # Minimum threshold for bound
    min_cdiffusion = 0.001  # Minimum threshold for constrained diffusion
    min_fdiffusion = 0.001  # Minimum threshold for fast diffusion



    # Toggle for including Gaussian fits in the combined plot
    include_gaussian_fits = False

    # Specify the input file path
    data = pd.read_csv(os.path.join(output_dir, ('gaps-and-fixes_decisions.csv' if configs['toggle']['use_gap_fixed'] else 'bound_decisions.csv')))
    bound_counts = data['Bound'].value_counts()
    counts_0 = bound_counts.get(0, 0)
    counts_1 = bound_counts.get(1, 0)
    counts_2 = bound_counts.get(2, 0)

    percentage_0 = (counts_0 / (counts_0 + counts_1 + counts_2)) * 100
    percentage_1 = (counts_1 / (counts_0 + counts_1 + counts_2)) * 100
    percentage_2 = (counts_2 / (counts_0 + counts_1 + counts_2)) * 100

    # Get the directory of the input file to save output files in the same directory

    # First pass: separate events based on "Bound"
    events = []

    # Iterate through each video
    for video in data['Video #'].unique():
        video_data = data[data['Video #'] == video]

        # Iterate through each cell
        for cell in video_data['Cell'].unique():
            cell_data = video_data[video_data['Cell'] == cell]

            # Iterate through each track
            for track in cell_data['Track'].unique():
                track_data = cell_data[cell_data['Track'] == track]
                track_length = len(track_data)

                # Create segments based on the number of frames per subdivision
                subdivs = [track_data.iloc[i:i + frames_per_subdivision] for i in
                           range(0, track_length, frames_per_subdivision)]

                for sub_id, sub in enumerate(subdivs):
                    if len(sub) < frames_per_subdivision:
                        continue  # Skip segments smaller than the specified number of frames

                    initial_frame = None
                    last_frame = None
                    current_bound = sub.iloc[0]['Bound']
                    positions = []
                    event_number = 1

                    # Iterate through each frame in the subdivision
                    for i in range(len(sub)):
                        frame_data = sub.iloc[i]

                        # Record initial frame
                        if initial_frame is None:
                            initial_frame = frame_data['Frame']

                        # Add all frames to positions, including -1
                        positions.append((frame_data['x'], frame_data['y']))

                        # Check if the bound value changes or if it's the last frame
                        if i == len(sub) - 1 or sub.iloc[i + 1]['Bound'] != current_bound:
                            # Record last frame
                            last_frame = frame_data['Frame']

                            # Save the event data
                            events.append({
                                "Video #": video,
                                "Video Name": frame_data['Video Name'],
                                "Cell": cell,
                                "Track": f"{track}_{sub_id + 1}",
                                "Event": event_number,
                                "Bound": current_bound,
                                "Initial Frame": initial_frame,
                                "Last Frame": last_frame,
                                "Positions": positions
                            })

                            # Reset for the next event
                            positions = []
                            if i < len(sub) - 1:
                                current_bound = sub.iloc[i + 1]['Bound']
                                event_number += 1
                                initial_frame = sub.iloc[i + 1]['Frame']

    # Second pass: calculate distances and MSDs
    final_results_with_correct_frames = []

    for event in events:
        positions = event["Positions"]
        df = pd.DataFrame(positions, columns=['x_coor', 'y_coor'])
        df['x_coor_um'] = df['x_coor'] * pixel_to_micron
        df['y_coor_um'] = df['y_coor'] * pixel_to_micron
        if len(df) >= 4:
            msd = calculate_MSD(df)
            D, alpha = fit_MSD(msd, time_interval)
        else:
            D, alpha = np.nan, np.nan
        frame_count = event["Last Frame"] - event["Initial Frame"] + 1  # Include gaps in the total count
        final_results_with_correct_frames.append({
            "Video #": event["Video #"],
            "Video Name": event["Video Name"],
            "Cell": event["Cell"],
            "Track": event["Track"],
            "Event": event["Event"],
            "Diffusion Coefficient": D,
            "Alpha": alpha,
            "Count": frame_count,
            "Bound": event["Bound"],
            "Initial Frame": event["Initial Frame"],
            "Last Frame": event["Last Frame"]
        })

    # Create a dataframe for the final results
    final_results_with_correct_frames_df = pd.DataFrame(final_results_with_correct_frames)

    # Sort the final results by Video #, Cell, Track, and Event
    sorted_final_results_with_correct_frames_df = final_results_with_correct_frames_df.sort_values(
        by=["Video #", "Cell", "Track", "Event"])

    # Save the results to a CSV file
    output_csv_path = os.path.join(output_dir, 'Diffusion_Coefficient_Calculation.csv')
    sorted_final_results_with_correct_frames_df.to_csv(output_csv_path, index=False)

    # Load the data from the provided CSV file
    data = pd.read_csv(output_csv_path)

    # Replace infinite values with NaN and then drop them
    data.replace([np.inf, -np.inf], np.nan, inplace=True)
    data.dropna(subset=['Diffusion Coefficient'], inplace=True)

    # Apply minimum thresholds and filter out rows based on the diffusion coefficient cutoffs for each bound state
    filtered_data_plot = data[
        ((data['Bound'] == 0) & (data['Diffusion Coefficient'] <= cutoff_fdiffusion) & (
                data['Diffusion Coefficient'] >= min_fdiffusion)) |
        ((data['Bound'] == 1) & (data['Diffusion Coefficient'] <= cutoff_cdiffusion) & (
                data['Diffusion Coefficient'] >= min_cdiffusion)) |
        ((data['Bound'] == 2) & (data['Diffusion Coefficient'] <= cutoff_bound) & (
                data['Diffusion Coefficient'] >= min_bound))
        ]

    # Filter data for plotting and separate based on the Bound
    filtered_bound_0_plot = filtered_data_plot[filtered_data_plot['Bound'] == 0]
    filtered_bound_1_plot = filtered_data_plot[filtered_data_plot['Bound'] == 1]
    filtered_bound_2_plot = filtered_data_plot[filtered_data_plot['Bound'] == 2]

    # Separate original data based on the Bound for proportion calculations
    filtered_bound_0 = data[(data['Bound'] == 0) & (data['Diffusion Coefficient'] >= min_fdiffusion)]
    filtered_bound_1 = data[(data['Bound'] == 1) & (data['Diffusion Coefficient'] >= min_cdiffusion)]
    filtered_bound_2 = data[(data['Bound'] == 2) & (data['Diffusion Coefficient'] >= min_bound)]

    # Calculate the total number of entries for proportion calculations
    total_entries = data['Count'].sum()

    # Calculate the frequency and percentage for each Bound including all data
    count_0 = filtered_bound_0['Count'].sum()
    count_1 = filtered_bound_1['Count'].sum()
    count_2 = filtered_bound_2['Count'].sum()

    percentage_0_filt = (count_0 / total_entries) * 100
    percentage_1_filt = (count_1 / total_entries) * 100
    percentage_2_filt = (count_2 / total_entries) * 100

    # Filter out rows where the Count is lower than 4 for plotting
    filtered_bound_0_plot = filtered_bound_0_plot[filtered_bound_0_plot['Count'] >= 4]
    filtered_bound_1_plot = filtered_bound_1_plot[filtered_bound_1_plot['Count'] >= 4]
    filtered_bound_2_plot = filtered_bound_2_plot[filtered_bound_2_plot['Count'] >= 4]

    # Create a list to hold the replicated values for plotting
    replicated_data = []
    replicated_bounds = []

    # Replicate Diffusion Coefficient values according to Count for each Bound
    for _, row in filtered_bound_0_plot.iterrows():
        replicated_data.extend([row['Diffusion Coefficient']] * int(row['Count']))
        replicated_bounds.extend([0] * int(row['Count']))

    for _, row in filtered_bound_1_plot.iterrows():
        replicated_data.extend([row['Diffusion Coefficient']] * int(row['Count']))
        replicated_bounds.extend([1] * int(row['Count']))

    for _, row in filtered_bound_2_plot.iterrows():
        replicated_data.extend([row['Diffusion Coefficient']] * int(row['Count']))
        replicated_bounds.extend([2] * int(row['Count']))

    # Create a DataFrame for plotting with the weights for histogram
    plot_data = pd.DataFrame({
        'Diffusion Coefficient': replicated_data,
        'Bound': replicated_bounds
    })

    # Define the number of bins for consistency in log scale
    base_num_bins = 60  # Base number of bins for the scaling factor of 1
    num_bins = int(base_num_bins * bin_scaling_factor)
    bins = np.logspace(np.log10(0.0001), np.log10(10), num=num_bins)  # Adjust bins based on scaling factor

    # Create a PDF file for saving the plots
    pdf_output_path = os.path.join(output_dir, 'Diffusion_Coefficient_Plots.pdf')
    with PdfPages(pdf_output_path) as pdf:
        # Create individual histograms for each Bound with Gaussian fit
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        mean_bound, std_bound, scale_bound = plot_gaussian_fit(
            plot_data[plot_data['Bound'] == 2]['Diffusion Coefficient'],
            '#66cc33', 'Bound', axes[0], bins)
        mean_cdiffusion, std_cdiffusion, scale_cdiffusion = plot_gaussian_fit(
            plot_data[plot_data['Bound'] == 1]['Diffusion Coefficient'], '#3366cc', 'C.Diffusion', axes[1], bins)
        mean_fdiffusion, std_fdiffusion, scale_fdiffusion = plot_gaussian_fit(
            plot_data[plot_data['Bound'] == 0]['Diffusion Coefficient'], '#cc3366', 'F.Diffusion', axes[2], bins)

        plt.tight_layout()
        pdf.savefig()  # Save the current figure into the PDF
        plt.close()

        # Create combined histogram with all Gaussian fits
        plt.figure(figsize=(15, 5))

        # Scale the histograms
        plt.hist(plot_data[plot_data['Bound'] == 0]['Diffusion Coefficient'], bins=bins, color='#cc3366', alpha=0.75,
                 label=f'Fast Diffusion ({percentage_0:.2f}%)', weights=np.full(len(plot_data[plot_data['Bound'] == 0]),
                                                                                percentage_0 / 100 / len(
                                                                                    plot_data[
                                                                                        plot_data['Bound'] == 0])))
        plt.hist(plot_data[plot_data['Bound'] == 1]['Diffusion Coefficient'], bins=bins, color='#3366cc', alpha=0.75,
                 label=f'Constrained Diffusion ({percentage_1:.2f}%)',
                 weights=np.full(len(plot_data[plot_data['Bound'] == 1]),
                                 percentage_1 / 100 / len(plot_data[plot_data['Bound'] == 1])))
        plt.hist(plot_data[plot_data['Bound'] == 2]['Diffusion Coefficient'], bins=bins, color='#66cc33', alpha=0.75,
                 label=f'Bound ({percentage_2:.2f}%)', weights=np.full(len(plot_data[plot_data['Bound'] == 2]),
                                                                       percentage_2 / 100 / len(
                                                                           plot_data[plot_data['Bound'] == 2])))

        plt.xscale('log')
        plt.xlim(0.0001, 50)  # Set limits from 0.0001 to 50 on the log scale
        plt.xticks([0.0001, 0.001, 0.01, 0.1, 1, 10], ['0.0001', '0.001', '0.01', '0.1', '1', '10'])  # Custom ticks
        plt.ylim(0, 0.07)

        # Plot Gaussian fits for each Bound if the toggle is on
        if include_gaussian_fits:
            x = np.linspace(np.log10(0.0001), np.log10(50), 100)
            if mean_bound is not None and std_bound is not None:
                p_bound = norm.pdf(x, mean_bound, std_bound)
                plt.plot(10 ** x, p_bound * scale_bound * (percentage_2 / 100), 'k', linewidth=2,
                         label='Bound Gaussian Fit')
            if mean_cdiffusion is not None and std_cdiffusion is not None:
                p_cdiffusion = norm.pdf(x, mean_cdiffusion, std_cdiffusion)
                plt.plot(10 ** x, p_cdiffusion * scale_cdiffusion * (percentage_1 / 100), 'k', linewidth=2,
                         label='C.Diffusion Gaussian Fit')
            if mean_fdiffusion is not None and std_fdiffusion is not None:
                p_fdiffusion = norm.pdf(x, mean_fdiffusion, std_fdiffusion)
                plt.plot(10 ** x, p_fdiffusion * scale_fdiffusion * (percentage_0 / 100), 'k', linewidth=2,
                         label='F.Diffusion Gaussian Fit')

        plt.title('Distribution of Apparent Diffusion Coefficient for All Behaviours with Gaussian Fits')
        plt.xlabel('Diffusion Coefficient (µm²/s)')
        plt.ylabel('Frequency')
        plt.legend(loc='upper left', bbox_to_anchor=(0, 1))
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y * 100:.0f}%'))
        pdf.savefig()  # Save the current figure into the PDF
        plt.close()

    print("Plots saved to PDF:", pdf_output_path)
    print("Results saved to CSV:", output_csv_path)
    return


# Function to calculate MSD
def calculate_MSD(df):
    N = len(df)
    msd = []
    for n in range(1, N):
        displacements = [(df['x_coor_um'].iloc[i + n] - df['x_coor_um'].iloc[i]) ** 2 +
                         (df['y_coor_um'].iloc[i + n] - df['y_coor_um'].iloc[i]) ** 2 for i in range(N - n)]
        msd.append(np.mean(displacements))
    return msd


# Function to fit MSD to extract the diffusion coefficient and alpha
def fit_MSD(msd, frame_interval):
    def msd_func(t, D, alpha):
        return 4 * D * (t ** alpha)

    times = np.arange(1, len(msd) + 1) * frame_interval
    popt, _ = curve_fit(msd_func, times, msd, bounds=(0, [np.inf, 2]))
    return popt


# Function to filter outliers using IQR method
def filter_outliers(data):
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return data[(data >= lower_bound) & (data <= upper_bound)]


# Function to plot Gaussian fit and convert y-axis to percentage
def plot_gaussian_fit(data, color, label, ax, bins):
    sns.histplot(data, kde=False, color=color, bins=bins, element='bars', stat='probability', ax=ax)
    ax.set_xscale('log')
    ax.set_xlim(0.0001, 50)  # Set limits from 0.0001 to 50 on the log scale
    ax.set_xticks([0.0001, 0.001, 0.01, 0.1, 1, 10])
    ax.set_xticklabels(['0.0001', '0.001', '0.01', '0.1', '1', '10'])
    ax.set_ylim([0, .15])

    # Convert y-axis to percentage
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y * 100:.0f}%'))

    # Fit Gaussian
    if len(data) > 0:
        log_data = np.log10(data)
        log_data = filter_outliers(log_data)  # Filter outliers
        mean, std = norm.fit(log_data)
        xmin, xmax = ax.get_xlim()
        x = np.linspace(np.log10(xmin), np.log10(xmax), 100)
        p = norm.pdf(x, mean, std)

        # Calculate the height of the peak bin
        hist, bin_edges = np.histogram(data, bins=bins, density=False)
        peak_bin_index = np.argmax(hist)
        peak_x = (bin_edges[peak_bin_index] + bin_edges[peak_bin_index + 1]) / 2
        peak_y = hist[peak_bin_index] / len(data)  # Normalize to the total count

        # Find the Gaussian value at peak_x
        peak_gaussian_y = norm.pdf(np.log10(peak_x), mean, std)

        # Scale factor to match the histogram's peak y value
        scale_factor = peak_y / peak_gaussian_y
        optimized_scale_factor, _ = curve_fit(lambda x, scale: norm.pdf(x, mean, std) * scale,
                                              np.log10(bins[:-1] + np.diff(bins) / 2), hist / len(data),
                                              p0=[scale_factor])

        ax.plot(10 ** x, p * optimized_scale_factor[0], 'k', linewidth=2)  # Scale the Gaussian fit
        ax.legend([f'{label} (peak D={10 ** mean:.6f} µm²/s)'])

        return mean, std, optimized_scale_factor[0]  # Return the Gaussian parameters
    else:
        ax.legend([f'{label} (no data available)'])
        return None, None, None

    ax.set_title(f'Distribution of Apparent Diffusion Coefficient for {label}')
    ax.set_xlabel('Diffusion Coefficient (µm²/s)')
    ax.set_ylabel('Frequency')


if __name__ == '__main__':
    start_time = time.time()
    parser = argparse.ArgumentParser(
        prog='rebind_MSD',
        description='#',
        epilog='Prior scripts: bound_classification.py, gaps_and_fixes.py')
    parser.add_argument('-c', '--config', default=None, type=str)
    args = parser.parse_args()
    print('Running as __main__ \n\t-> config_path: ' + str(args.config))
    main(args.config)
    print("--- %s seconds ---" % (time.time() - start_time))