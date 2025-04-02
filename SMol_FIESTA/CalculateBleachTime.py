import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import expon
import warnings

warnings.filterwarnings("ignore")

#%% PARAMETERS

# Path to the input CSV file (change as needed)
input_csv = r'D:\Microscopy\SMSNEW\_ANIMAL_CELL_DATA\Rif1\pablo\im\ImagesNoMetadata\AnalysisRebindCBC_start0_Quality3n5_noConFilt\SF\gaps-and-fixes_decisions.csv'
# Output CSV for the grouped data (unchanged)
output_csv = 'processed_video_data_with_frames.csv'
# Output PDF for the plot
output_pdf = r'D:\Microscopy\SMSNEW\_ANIMAL_CELL_DATA\Rif1\pablo\im\ImagesNoMetadata\AnalysisRebindCBC_start0_Quality3n5_noConFilt\SF\b_plot.pdf'

# Scaling factor for length (e.g. if you want to divide by 100 as in your original script)
scaling_factor = 1/100.0

# Bootstrapping parameters
n_bootstrap = 1000  # number of bootstrap samples

#%% LOAD AND GROUP DATA

print("Loading data...")
data = pd.read_csv(input_csv)

# Group the data by Video #, Video Name, Cell, Track
grouped = data.groupby(['Video #', 'Video Name', 'Cell', 'Track']).agg({'Frame': ['min', 'max']}).reset_index()
grouped.columns = ['Video #', 'Video Name', 'Cell', 'Track', 'Min_Frame', 'Max_Frame']

# Compute track Length as max - min, then scale
grouped['Length'] = (grouped['Max_Frame'] - grouped['Min_Frame']) * scaling_factor

# Save the grouped data (this output remains the same)
grouped[['Video #', 'Video Name', 'Cell', 'Track', 'Length', 'Min_Frame', 'Max_Frame']].to_csv(output_csv, index=False)
print("Processed grouped data saved to", output_csv)

#%% EXPONENTIAL MODEL & MLE ESTIMATION WITH LEFT-TRUNCATION

# We assume that track lengths follow an exponential decay but we can only observe t >= L_min.
# Define L_min as the minimum observed length.
lengths = grouped['Length'].values
L_min = lengths.min()

# Only use lengths above the minimum (they are all >= L_min by definition)
# The MLE estimator for lambda is: λ_hat = n / sum(t_i - L_min)
n = len(lengths)
lambda_hat = n / np.sum(lengths - L_min)
mean_estimated = L_min + 1.0/lambda_hat

print("\nRobust MLE estimation (using left-truncated exponential):")
print("Number of tracks (n):", n)
print("Minimum observed length (L_min):", L_min)
print("Estimated lambda: {:.4f}".format(lambda_hat))
print("Estimated mean track length (extrapolated): {:.4f}".format(mean_estimated))

#%% BOOTSTRAPPING FOR UNCERTAINTY

def mle_mean_estimator(sample, L_min):
    n_sample = len(sample)
    lam_hat = n_sample / np.sum(sample - L_min)
    return L_min + 1.0/lam_hat

bootstrap_means = []
rng = np.random.default_rng()  # use numpy random generator
for i in range(n_bootstrap):
    sample = rng.choice(lengths, size=n, replace=True)
    bootstrap_means.append(mle_mean_estimator(sample, L_min))
bootstrap_means = np.array(bootstrap_means)
mean_bootstrap = bootstrap_means.mean()
se_bootstrap = bootstrap_means.std()
ci_lower = np.percentile(bootstrap_means, 2.5)
ci_upper = np.percentile(bootstrap_means, 97.5)

print("\nBootstrap results ({} iterations):".format(n_bootstrap))
print("Bootstrap mean: {:.4f}".format(mean_bootstrap))
print("Standard error: {:.4f}".format(se_bootstrap))
print("95% CI: [{:.4f}, {:.4f}]".format(ci_lower, ci_upper))

#%% PLOTTING

# Define the theoretical PDF for the truncated exponential:
# f(t) = lambda * exp(-lambda*(t - L_min))  for t >= L_min.
def truncated_exponential_pdf(t, lam, L_min):
    return lam * np.exp(-lam*(t - L_min))

# Create the histogram
fig, ax = plt.subplots(figsize=(8, 6))
bins = np.linspace(L_min, lengths.max(), 50)
hist_counts, bin_edges, _ = ax.hist(lengths, bins=bins, density=True, alpha=0.5, label='Observed Data')

# Create x-values for plotting the fitted PDF
x_vals = np.linspace(L_min, lengths.max(), 200)
fitted_pdf = truncated_exponential_pdf(x_vals, lambda_hat, L_min)

ax.plot(x_vals, fitted_pdf, 'r-', linewidth=2, label='Fitted Exponential')
ax.axvline(mean_estimated, color='green', linestyle='--', linewidth=2,
           label=f'Extrapolated Mean: {mean_estimated:.2f}')

# Annotate the plot with the estimated parameters
textstr = '\n'.join((
    f"n = {n}",
    f"L_min = {L_min:.2f}",
    f"λ̂ = {lambda_hat:.4f}",
    f"Mean = {mean_estimated:.2f}",
    f"Bootstrap SE = {se_bootstrap:.2f}",
    f"95% CI = [{ci_lower:.2f}, {ci_upper:.2f}]"
))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

ax.set_xlabel('Track Length (scaled)')
ax.set_ylabel('Probability Density')
ax.set_title('Histogram of Track Lengths with Fitted Exponential Decay')
ax.legend()

plt.tight_layout()
plt.savefig(output_pdf, format='pdf', bbox_inches='tight')
plt.close()
print("\nPDF plot saved as", output_pdf)
