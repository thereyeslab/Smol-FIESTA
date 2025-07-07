import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import expon
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
from sklearn.linear_model import LogisticRegression
import tomllib
import argparse
import time
import matplotlib.ticker as ticker

import matplotlib as mpl

mpl.rcParams.update({
    # use Arial for all text
    'font.family':     'sans-serif',
    'font.sans-serif': ['Arial'],
    # embed fonts as Type 42 (TrueType) so Illustrator can edit them
    'pdf.fonttype':    42,
    'ps.fonttype':     42,
    # optionally turn off the grid globally
    'axes.grid':       False,
})


# Error handling for threadpoolctl issue
os.environ["OMP_NUM_THREADS"] = "1"
def train_logistic_on_synthetic(lambda1, lambda2, a=4.0, n=20000):
    """Simulate and train a logistic model to classify Exp1 vs Exp2."""
    x1 = np.random.exponential(1 / lambda1, size=n) + a
    x2 = np.random.exponential(1 / lambda2, size=n) + a
    X  = np.concatenate([x1, x2])
    y  = np.concatenate([np.ones(n), np.zeros(n)])

    x_adj = X - a
    S1    = np.exp(-lambda1 * a)
    S2    = np.exp(-lambda2 * a)
    f1    = lambda1 * np.exp(-lambda1 * x_adj) / S1
    f2    = lambda2 * np.exp(-lambda2 * x_adj) / S2

    eps = 1e-300
    llr = np.log(np.clip(f1, eps, None) / np.clip(f2, eps, None)).reshape(-1, 1)

    model = LogisticRegression(solver='lbfgs')
    model.fit(llr, y)
    return model

def neg_log_likelihood_single(params, data, a=4.0):
    lambda1 = params[0]
    x_adj = data - a
    if np.any(x_adj < 0):
        return np.inf  # invalid under truncation
    S = np.exp(-lambda1 * a)
    likelihoods = lambda1 * np.exp(-lambda1 * x_adj) / S
    neg_log_lik = -np.sum(np.log(likelihoods))
    return neg_log_lik


def neg_log_likelihood_with_bleaching(params, data, lambda_bleach):
    lambda1, lambda2, p1 = params
    p2 = 1 - p1

    # Combined decay rates including photobleaching
    lambda1_total = lambda1
    lambda2_total = lambda2

    # Adjusted likelihoods for the mixture of two exponentials with photobleaching
    likelihoods = (p1 * expon.pdf(data, scale=1 / lambda1_total) +
                   p2 * expon.pdf(data, scale=1 / lambda2_total))

    neg_log_lik = -np.sum(np.log(likelihoods))
    return neg_log_lik
def neg_log_likelihood_trunc(params, data, a):
    lambda1, lambda2, p1 = params
    p2 = 1 - p1

    x_adj = data - a
    if np.any(x_adj < 0):
        return np.inf

    # Do NOT divide by S1 or S2 — just use truncated form directly
    f1 = p1 * lambda1 * np.exp(-lambda1 * x_adj)
    f2 = p2 * lambda2 * np.exp(-lambda2 * x_adj)

    likelihoods = f1 + f2
    likelihoods = np.clip(likelihoods, 1e-300, None)  # prevent log(0)

    return -np.sum(np.log(likelihoods))


def get_initial_parameters(data):
    data_reshaped = data.reshape(-1, 1)
    kmeans = KMeans(n_clusters=2, n_init=10, random_state=0).fit(data_reshaped)
    labels = kmeans.labels_
    cluster_centers = kmeans.cluster_centers_.flatten()
    lambda1_initial = 1 / cluster_centers[0]
    lambda2_initial = 1 / cluster_centers[1]
    p1_initial = np.sum(labels == 0) / len(data)
    return lambda1_initial, lambda2_initial, p1_initial

def run_single_exponential_model(data, a=4.0):
    x_adj = data - a
    initial_lambda = 1 / np.mean(x_adj)
    result = minimize(neg_log_likelihood_single, [initial_lambda], args=(data, a), bounds=[(1e-6, None)])
    lambda1 = result.x[0]
    neg_log_lik = neg_log_likelihood_single([lambda1], data, a)
    n = len(data)
    bic = 2 * neg_log_lik + np.log(n) * 1
    mean_time_single = 1 / lambda1 + a
    return lambda1, mean_time_single, bic


def run_mixture_model(data_path, event_type=None, max_iterations=100, tolerance=0.05, frameRate=1):
    # load and subset
    df = pd.read_csv(data_path)
    if event_type:
        raw = df.loc[df['type'] == event_type, 'time'].values
    else:
        raw = df['time'].values
    data = raw / frameRate

    # fixed left-truncation at 4.0
    a = np.min(data)

    # EM-like fit of truncated double exponential
    for _ in range(max_iterations):
        l1_init, l2_init, p1_init = get_initial_parameters(data)
        result = minimize(
            neg_log_likelihood_trunc,
            [l1_init, l2_init, p1_init],
            args=(data, a),
            bounds=[(1e-6, None), (1e-6, None), (0, 1)]
        )
        l1, l2, p1 = result.x
        p2 = 1 - p1
        # ensure λ1 <-> larger mean
        if l1 < l2:
            l1, l2 = l2, l1
            p1, p2 = p2, p1

    # compute BIC for double exp
    neg_ll = neg_log_likelihood_trunc([l1, l2, p1], data, a)
    n = len(data)
    bic_double = 2 * neg_ll + np.log(n) * 2  # two decay rates + one weight

    # corrected vs noncorrected means
    mean1_nc = 1 / l1
    mean2_nc = 1 / l2
    label = 'Bound time' if event_type == 'Bound' else 'Diffusion time'
    # assemble results
    results = { label:"",
        "Estimated Mean Time 1 (non-corrected)": mean1_nc + a,
        "Estimated Decay Rate 1 (non-corrected)": l1,
        "Estimated Mean Time 2 (non-corrected)": mean2_nc + a,
        "Estimated Decay Rate 2 (non-corrected)": l2,
        "Fraction of Exp Component 1":           p1,
        "Fraction of Exp Component 2":           p2,
        "Weight Fraction 1":                     (mean1_nc * p1) / (mean1_nc * p1 + mean2_nc * p2),
        "Weight Fraction 2":                     (mean2_nc * p2) / (mean1_nc * p1 + mean2_nc * p2),
        "BIC Double Exponential":                bic_double
    }
    return results, data, l1, l2, p1, p2

def plot_results(data, λ1, λ2, p1, p2, λ1_single=None, mean_time_single=None):
    """
    Big-text, no-grid histogram with single or double exponential fit.
    """
    a       = data.min()
    bins    = np.linspace(0, data.max(), 50)
    hist, _ = np.histogram(data, bins=bins, density=True)
    centers = 0.5 * (bins[1:] + bins[:-1])

    plt.figure()
    plt.grid(False)

    # raw data
    plt.hist(data, bins=bins, density=True,
             color='lightgray', edgecolor='black', alpha=0.6)

    if λ1_single is not None:
        # single-exp
        pdf = expon.pdf(centers, scale=1/λ1_single)
        plt.plot(centers, pdf, 'g--', linewidth=3,
                 label=f'Single exp\n(mean {(1/λ1_single)+a:.2f}s)')
    else:
        # double-exp filled bands + total fit
        pdf1 = p1        * expon.pdf(centers, scale=1/λ1)
        pdf2 = (1.0 - p1) * expon.pdf(centers, scale=1/λ2)
        plt.fill_between(centers, 0,       pdf2,       color='C1', alpha=0.4,
                         label=f'2nd exp (mean {(1/λ2)+a:.2f}s)')
        plt.fill_between(centers, pdf2,    pdf1 + pdf2, color='C0', alpha=0.4,
                         label=f'1st exp (mean {(1/λ1)+a:.2f}s)')
        plt.plot(centers, pdf1 + pdf2, 'k-', linewidth=3, label='Total fit')

    # force x-axis to start at zero
    ax = plt.gca()
    ax.set_xlim(left=0)

    # big, sparse ticks
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.tick_params(axis='both', which='major',
                   labelsize=20, length=10, width=2)

    # labels & legend
    plt.xlabel('time (seconds)', fontsize=24)
    plt.ylabel('PDF',              fontsize=24)
    plt.title('Histogram with Exponential Fits', fontsize=28)
    plt.legend(fontsize=20, frameon=False)
    plt.tight_layout()


def plot_results(data, λ1, λ2, p1, p2, λ1_single=None, mean_time_single=None):
    """
    Big‐text, no‐grid histogram with single or double exponential fit.
    """
    a       = data.min()
    bins    = np.linspace(0, data.max(), 50)
    hist, _ = np.histogram(data, bins=bins, density=True)
    centers = 0.5 * (bins[1:] + bins[:-1])

    plt.figure(figsize=(8,6))
    plt.grid(False)

    # gray bars
    plt.hist(data, bins=bins, density=True,
             color='#858282', edgecolor='black', alpha=0.6)

    if λ1_single is not None:
        # single‐exp fit
        pdf = expon.pdf(centers, scale=1/λ1_single)
        plt.plot(centers, pdf, 'g--', linewidth=3,
                 label=f'Single exp\n(mean {(1/λ1_single)+a:.2f}s)')
    else:
        # double‐exp filled bands
        pdf1 = p1        * expon.pdf(centers, scale=1/λ1)
        pdf2 = (1.0 - p1) * expon.pdf(centers, scale=1/λ2)
        # 2nd exp band in blue
        plt.fill_between(centers, 0, pdf2,
                         color='#1f77b4', alpha=0.4,
                         label=f'2nd exp (mean {(1/λ2)+a:.2f}s)')
        # 1st exp band in red
        plt.fill_between(centers, pdf2, pdf1+pdf2,
                         color='r', alpha=0.4,
                         label=f'1st exp (mean {(1/λ1)+a:.2f}s)')
        plt.plot(centers, pdf1+pdf2, 'k-', linewidth=3,
                 label='Total fit')

    # force x axis to start at zero
    ax = plt.gca()
    ax.set_xlim(left=0)

    # big, sparse ticks
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.tick_params(axis='both', which='major',
                   labelsize=24, length=10, width=2)

    # labels & legend
    plt.xlabel('time (seconds)', fontsize=24)
    plt.ylabel('PDF',              fontsize=26)
    plt.title('Histogram with Exponential Fits', fontsize=28)
    plt.legend(fontsize=24, frameon=False)
    plt.tight_layout()


def plot_results_with_trendlines(data, λ1, λ2, p1, p2, λ1_single=None, mean_time_single=None):
    """
    Big‐text, no‐grid histogram with dashed trendlines only.
    """
    a       = data.min()
    bins    = np.linspace(0, data.max(), 50)
    hist, _ = np.histogram(data, bins=bins, density=True)
    centers = 0.5 * (bins[1:] + bins[:-1])

    plt.figure(figsize=(8,6))
    plt.grid(False)

    # gray bars
    plt.hist(data, bins=bins, density=True,
             color='#858282', edgecolor='black', alpha=0.6)

    # dashed trendlines
    pdf1 = p1        * expon.pdf(centers, scale=1/λ1)
    pdf2 = (1.0 - p1) * expon.pdf(centers, scale=1/λ2)
    plt.plot(centers, pdf1, 'r--', linewidth=3,
             label=f'1st exp\n(mean {(1/λ1)+a:.2f}s)')
    plt.plot(centers, pdf2, color='#1f77b4', linestyle='--', linewidth=3,
             label=f'2nd exp\n(mean {(1/λ2)+a:.2f}s)')

    # force x axis to start at zero
    ax = plt.gca()
    ax.set_xlim(left=0)

    # big, sparse ticks
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.tick_params(axis='both', which='major',
                   labelsize=24, length=10, width=2)

    # labels & legend
    plt.xlabel('time (seconds)', fontsize=24)
    plt.ylabel('PDF',              fontsize=26)
    plt.title('Histogram with Double Exponential Trendlines', fontsize=28)
    plt.legend(fontsize=24, frameon=False)
    plt.tight_layout()



def main(config_path:str = None):
    if not config_path:
        __location__ = os.path.realpath(
            os.path.join(os.getcwd(), os.path.dirname(__file__)))
        config_path = os.path.join(__location__, 'script-config.toml')
    with open(config_path, 'rb') as config_file:
        configs = tomllib.load(config_file)

    csv_path = configs['path']['csv_path']
    output_path = str(os.path.join(csv_path, configs['path']['output_folder_name']))
    data_path = str(os.path.join(output_path, 'rebind-Events.csv'))

    # Path, measurement and frame rate
    interval    = configs['MSD-analysis']['time_interval_ms']
    frameRate   = 1000 / interval

    # Toggles
    process_bound  = configs['two-exponential-analysis']['Test_Bound']
    process_search = configs['two-exponential-analysis']['Test_Diffusion']
    output_folder = os.path.dirname(data_path)
    pdf_path      = os.path.join(output_folder, '2Exponential_test.pdf')

    # load once and prepare merged posterior column
    raw_df = pd.read_csv(data_path)
    raw_df['Posterior_Prob_Exp1'] = "NA"
    with PdfPages(pdf_path) as pdf:
        for event_type, test_double in [('Bound', process_bound), ('SearchTime', process_search)]:
            data_vals = raw_df.loc[raw_df['type'] == event_type, 'time'].values / frameRate
            a = np.min(data_vals)
            lambda1_s, mean_s, bic_s = run_single_exponential_model(data_vals, a=a)
            xlabel = 'Bound time' if event_type == 'Bound' else 'Diffusion time'

            if test_double:
                results, data_vals, lambda1, lambda2, p1, p2 = run_mixture_model(
                    data_path=data_path,
                    event_type=event_type,
                    max_iterations=100,
                    tolerance=0.05,
                    frameRate=frameRate
                )
                bic_d = results["BIC Double Exponential"]

                if bic_s < bic_d:
                    plt.figure(figsize=(8, 4))
                    plt.text(
                        0.1, 0.1,
                        f"{event_type} single-exp\nMean={mean_s:.3f}, λ={lambda1_s:.3f}, BIC single={bic_s:.1f}",
                        fontsize=12
                    )
                    plt.axis('off');
                    pdf.savefig();
                    plt.close()

                    plot_results(data_vals, None, None, None, None, lambda1_s, mean_s)
                    plt.xlabel(xlabel);
                    pdf.savefig();
                    plt.close()

                    raw_df.loc[raw_df['type'] == event_type, 'Posterior_Prob_Exp1'] = 1.0

                else:
                    plt.figure(figsize=(8, 4))
                    plt.text(
                        0.1, 0.05,
                        '\n'.join(f"{k}: {v}" for k, v in results.items()),
                        fontsize=12
                    )
                    plt.axis('off');
                    pdf.savefig();
                    plt.close()

                    plot_results(data_vals, lambda1, lambda2, p1, p2)
                    plt.xlabel(xlabel);
                    pdf.savefig();
                    plt.close()

                    plot_results_with_trendlines(data_vals, lambda1, lambda2, p1, p2)
                    plt.xlabel(xlabel);
                    pdf.savefig();
                    plt.close()

                    logistic_model = train_logistic_on_synthetic(lambda1, lambda2, a=a)
                    x_adj = data_vals - a
                    S1 = np.exp(-lambda1 * a)
                    S2 = np.exp(-lambda2 * a)
                    f1 = lambda1 * np.exp(-lambda1 * x_adj) / S1
                    f2 = lambda2 * np.exp(-lambda2 * x_adj) / S2
                    llr = np.log(np.clip(f1, 1e-300, None) / np.clip(f2, 1e-300, None)).reshape(-1, 1)
                    post = logistic_model.predict_proba(llr)[:, 1]
                    raw_df.loc[raw_df['type'] == event_type, 'Posterior_Prob_Exp1'] = post

            else:
                plt.figure(figsize=(8, 4))
                plt.text(
                    0.1, 0.1,
                    f"{event_type} forced single-exp\nMean={mean_s:.3f}, λ={lambda1_s:.3f}, BIC={bic_s:.1f}",
                    fontsize=12
                )
                plt.axis('off');
                pdf.savefig();
                plt.close()

                plot_results(data_vals, None, None, None, None, lambda1_s, mean_s)
                plt.xlabel(xlabel);
                pdf.savefig();
                plt.close()

    # write out
    print(f"Results and plots saved to {pdf_path}")
    merged_csv = os.path.join(output_folder,"Intermediates", 'Exp1_probabilities.csv')
    raw_df.to_csv(merged_csv, index=False)
    print(f"Merged posterior CSV saved to {merged_csv}")

    # === Transition matrix calculation ===

    # 1) Define which event types to include
    event_types = ['Bound', 'SearchTime']
    present     = [t for t in event_types if t in raw_df['type'].unique()]

    # 2) Detect number of components per type
    exp_counts = {}
    for t in present:
        p_vals = raw_df.loc[raw_df['type'] == t, 'Posterior_Prob_Exp1'] \
                       .dropna() \
                       .astype(float)
        # if any posterior < 1 → two exps, else one
        exp_counts[t] = 2 if (p_vals < 0.999999).any() else 1

    # 3) Build state labels and index map
    state_labels = []
    for t in present:
        for comp in range(1, exp_counts[t] + 1):
            state_labels.append(f"{t}_Exp{comp}")
    state_idx = {s: i for i, s in enumerate(state_labels)}

    # 4) Initialize counts matrix
    n_states     = len(state_labels)
    trans_counts = np.zeros((n_states, n_states), dtype=float)

    # 5) Filter to only those types and group
    df2     = raw_df[raw_df['type'].isin(present)].copy()
    grouped = df2.groupby(['Video #', 'Cell', 'Track'])

    for _, grp in grouped:
        grp = grp.sort_values('StartFrame')
        prev_probs = None
        prev_type  = None

        for _, row in grp.iterrows():
            t    = row['type']
            dur  = row['time']
            p1   = float(row['Posterior_Prob_Exp1'])
            n    = exp_counts[t]
            # soft‐assignment vector
            probs = [p1, 1 - p1] if n == 2 else [1.0]

            # a) self‐transitions (diagonal)
            for comp, pr in enumerate(probs, start=1):
                i = state_idx[f"{t}_Exp{comp}"]
                trans_counts[i, i] += dur * pr

            # b) cross‐transitions from previous event
            if prev_type is not None:
                for comp_i, pi in enumerate(prev_probs, start=1):
                    i = state_idx[f"{prev_type}_Exp{comp_i}"]
                    for comp_j, pj in enumerate(probs, start=1):
                        j = state_idx[f"{t}_Exp{comp_j}"]
                        trans_counts[i, j] += pi * pj

            prev_type, prev_probs = t, probs

    # 6) Row‐normalize → probabilities
    row_sums    = trans_counts.sum(axis=1, keepdims=True)
    trans_probs = np.divide(
        trans_counts,
        row_sums,
        out=np.zeros_like(trans_counts),
        where=row_sums != 0
    )

    # 7) Save & display
    df_trans  = pd.DataFrame(trans_probs, index=state_labels, columns=state_labels)
    print(df_trans)
    out_path  = os.path.join(output_folder, 'transition_probabilities.csv')
    df_trans.to_csv(out_path, float_format='%.6f')
    print(f"Transition probabilities saved to {out_path}")


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
