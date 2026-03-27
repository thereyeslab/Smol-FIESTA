from __future__ import annotations

import os
import time
import logging
import argparse
import tomllib
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class BehaviorBootstrapResult:
    point_estimates: dict[str, float]
    bootstrap_means: dict[str, float]
    standard_errors: dict[str, float]
    ci_95: dict[str, tuple[float, float]]
    n_boot: int
    n_units: int
    unit_name: str
    bootstrap_estimates: pd.DataFrame


BEHAVIOR_ORDER = ["F. Diffusion", "C. Diffusion", "Bound"]


def logging_setup(path: str, script_name: str) -> str:
    log_dir = os.path.join(path, "logs")
    os.makedirs(log_dir, exist_ok=True)

    log_file = os.path.join(log_dir, f"LOG_{script_name}.txt")
    handlers = [logging.FileHandler(log_file)]
    logging.basicConfig(
        format="%(message)s",
        level=logging.INFO,
        handlers=handlers,
        force=True,
    )

    open(log_file, "w").close()
    return log_file


def print_log(*args, end: str = "\n") -> None:
    msg = " ".join(str(a) for a in args)
    print(msg, end=end)
    logging.info(msg + end)


def validate_spot_table(df: pd.DataFrame) -> None:
    required_cols = {"Video #", "Cell", "Track", "Frame", "Bound"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if df.empty:
        raise ValueError("Input dataframe is empty.")


def _count_label(x: pd.Series, label: int) -> int:
    return int((x == label).sum())


def prepare_track_behavior_summary(spots_df: pd.DataFrame) -> pd.DataFrame:
    validate_spot_table(spots_df)

    df = spots_df.copy()
    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]
    df = df.sort_values(["Video #", "Cell", "Track", "Frame"]).reset_index(drop=True)

    group_cols = ["Video #", "Cell", "Track"]

    track_df = (
        df.groupby(group_cols, as_index=False)
        .agg(
            total_frames=("Bound", "size"),
            n_fast_diffusion=("Bound", lambda x: _count_label(x, 0)),
            n_constrained_diffusion=("Bound", lambda x: _count_label(x, 1)),
            n_bound=("Bound", lambda x: _count_label(x, 2)),
        )
    )

    check_sum = (
        track_df["n_fast_diffusion"]
        + track_df["n_constrained_diffusion"]
        + track_df["n_bound"]
    )
    if not (check_sum == track_df["total_frames"]).all():
        raise ValueError(
            "Behavior counts do not sum to total_frames for one or more tracks."
        )

    return track_df


def prepare_cell_behavior_summary(track_df: pd.DataFrame) -> pd.DataFrame:
    required_cols = {
        "Video #",
        "Cell",
        "total_frames",
        "n_fast_diffusion",
        "n_constrained_diffusion",
        "n_bound",
    }
    missing = required_cols - set(track_df.columns)
    if missing:
        raise ValueError(
            f"Track summary missing required columns for cell aggregation: {sorted(missing)}"
        )

    cell_df = (
        track_df.groupby(["Video #", "Cell"], as_index=False)
        .agg(
            total_frames=("total_frames", "sum"),
            n_fast_diffusion=("n_fast_diffusion", "sum"),
            n_constrained_diffusion=("n_constrained_diffusion", "sum"),
            n_bound=("n_bound", "sum"),
        )
    )

    check_sum = (
        cell_df["n_fast_diffusion"]
        + cell_df["n_constrained_diffusion"]
        + cell_df["n_bound"]
    )
    if not (check_sum == cell_df["total_frames"]).all():
        raise ValueError(
            "Behavior counts do not sum to total_frames for one or more cells."
        )

    return cell_df


def compute_pooled_behavior_proportions(summary_df: pd.DataFrame) -> dict[str, float]:
    required_cols = {
        "total_frames",
        "n_fast_diffusion",
        "n_constrained_diffusion",
        "n_bound",
    }
    missing = required_cols - set(summary_df.columns)
    if missing:
        raise ValueError(
            f"Summary dataframe missing required columns: {sorted(missing)}"
        )

    total_frames = int(summary_df["total_frames"].sum())
    if total_frames == 0:
        return {name: 0.0 for name in BEHAVIOR_ORDER}

    fast_count = int(summary_df["n_fast_diffusion"].sum())
    constrained_count = int(summary_df["n_constrained_diffusion"].sum())
    bound_count = int(summary_df["n_bound"].sum())

    return {
        "F. Diffusion": fast_count / total_frames,
        "C. Diffusion": constrained_count / total_frames,
        "Bound": bound_count / total_frames,
    }


def _summarize_bootstrap_estimates(
    estimates_df: pd.DataFrame,
    point_estimates: dict[str, float],
    n_boot: int,
    n_units: int,
    unit_name: str,
) -> BehaviorBootstrapResult:
    bootstrap_means = {
        col: float(estimates_df[col].mean())
        for col in BEHAVIOR_ORDER
    }
    standard_errors = {
        col: float(estimates_df[col].std(ddof=1))
        for col in BEHAVIOR_ORDER
    }
    ci_95 = {
        col: (
            float(np.percentile(estimates_df[col], 2.5)),
            float(np.percentile(estimates_df[col], 97.5)),
        )
        for col in BEHAVIOR_ORDER
    }

    return BehaviorBootstrapResult(
        point_estimates=point_estimates,
        bootstrap_means=bootstrap_means,
        standard_errors=standard_errors,
        ci_95=ci_95,
        n_boot=n_boot,
        n_units=n_units,
        unit_name=unit_name,
        bootstrap_estimates=estimates_df,
    )


def bootstrap_behavior_proportions_by_track(
    spots_df: pd.DataFrame,
    n_boot: int = 1000,
    random_state: Optional[int] = 42,
) -> BehaviorBootstrapResult:
    track_df = prepare_track_behavior_summary(spots_df)

    if track_df.empty:
        raise ValueError("No tracks found after summarization.")

    rng = np.random.default_rng(random_state)
    n_tracks = len(track_df)

    point_estimates = compute_pooled_behavior_proportions(track_df)
    estimates = np.empty((n_boot, 3), dtype=float)

    for i in range(n_boot):
        sample_idx = rng.integers(0, n_tracks, size=n_tracks)
        sampled_tracks = track_df.iloc[sample_idx]
        props = compute_pooled_behavior_proportions(sampled_tracks)

        estimates[i, 0] = props["F. Diffusion"]
        estimates[i, 1] = props["C. Diffusion"]
        estimates[i, 2] = props["Bound"]

    estimates_df = pd.DataFrame(estimates, columns=BEHAVIOR_ORDER)

    return _summarize_bootstrap_estimates(
        estimates_df=estimates_df,
        point_estimates=point_estimates,
        n_boot=n_boot,
        n_units=n_tracks,
        unit_name="track",
    )


def bootstrap_behavior_proportions_by_cell(
    spots_df: pd.DataFrame,
    n_boot: int = 1000,
    random_state: Optional[int] = 42,
) -> BehaviorBootstrapResult:
    track_df = prepare_track_behavior_summary(spots_df)
    cell_df = prepare_cell_behavior_summary(track_df)   
    cell_df["bound_fraction"] = cell_df["n_bound"] / cell_df["total_frames"]
    print(cell_df["bound_fraction"].describe())
    import matplotlib.pyplot as plt
    plt.hist(cell_df["bound_fraction"], bins=10)
    plt.show()

    if cell_df.empty:
        raise ValueError("No cells found after summarization.")

    rng = np.random.default_rng(random_state)
    n_cells = len(cell_df)

    point_estimates = compute_pooled_behavior_proportions(cell_df)
    estimates = np.empty((n_boot, 3), dtype=float)

    for i in range(n_boot):
        sample_idx = rng.integers(0, n_cells, size=n_cells)
        sampled_cells = cell_df.iloc[sample_idx]
        props = compute_pooled_behavior_proportions(sampled_cells)

        estimates[i, 0] = props["F. Diffusion"]
        estimates[i, 1] = props["C. Diffusion"]
        estimates[i, 2] = props["Bound"]

    estimates_df = pd.DataFrame(estimates, columns=BEHAVIOR_ORDER)

    return _summarize_bootstrap_estimates(
        estimates_df=estimates_df,
        point_estimates=point_estimates,
        n_boot=n_boot,
        n_units=n_cells,
        unit_name="cell",
    )


def bootstrap_behavior_proportions_hierarchical(
    spots_df: pd.DataFrame,
    n_boot: int = 1000,
    random_state: Optional[int] = 42,
) -> BehaviorBootstrapResult:
    """
    Hierarchical bootstrap:
    1. sample cells with replacement
    2. within each sampled cell, sample its tracks with replacement
    3. recompute pooled frame-level proportions
    """
    track_df = prepare_track_behavior_summary(spots_df)

    if track_df.empty:
        raise ValueError("No tracks found after summarization.")

    rng = np.random.default_rng(random_state)

    grouped_cells = {
        (video, cell): group.reset_index(drop=True)
        for (video, cell), group in track_df.groupby(["Video #", "Cell"], sort=False)
    }

    cell_ids = list(grouped_cells.keys())
    n_cells = len(cell_ids)

    point_estimates = compute_pooled_behavior_proportions(track_df)
    estimates = np.empty((n_boot, 3), dtype=float)

    for i in range(n_boot):
        sampled_cell_ids = rng.choice(n_cells, size=n_cells, replace=True)

        sampled_tracks_per_replicate = []

        for cell_idx in sampled_cell_ids:
            cell_group = grouped_cells[cell_ids[cell_idx]]
            n_tracks_in_cell = len(cell_group)

            sampled_track_idx = rng.integers(
                0,
                n_tracks_in_cell,
                size=n_tracks_in_cell,
            )
            sampled_cell_tracks = cell_group.iloc[sampled_track_idx]
            sampled_tracks_per_replicate.append(sampled_cell_tracks)

        sampled_tracks = pd.concat(sampled_tracks_per_replicate, ignore_index=True)
        props = compute_pooled_behavior_proportions(sampled_tracks)

        estimates[i, 0] = props["F. Diffusion"]
        estimates[i, 1] = props["C. Diffusion"]
        estimates[i, 2] = props["Bound"]

    estimates_df = pd.DataFrame(estimates, columns=BEHAVIOR_ORDER)

    return _summarize_bootstrap_estimates(
        estimates_df=estimates_df,
        point_estimates=point_estimates,
        n_boot=n_boot,
        n_units=n_cells,
        unit_name="hierarchical_cell_track",
    )


def result_to_dataframe(result: BehaviorBootstrapResult) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Behavior": BEHAVIOR_ORDER,
            "PointEstimate": [result.point_estimates[k] for k in BEHAVIOR_ORDER],
            "BootstrapMean": [result.bootstrap_means[k] for k in BEHAVIOR_ORDER],
            "SE": [result.standard_errors[k] for k in BEHAVIOR_ORDER],
            "CI_Lower": [result.ci_95[k][0] for k in BEHAVIOR_ORDER],
            "CI_Upper": [result.ci_95[k][1] for k in BEHAVIOR_ORDER],
            "BootstrapUnit": result.unit_name,
            "N_Units": result.n_units,
            "N_Boot": result.n_boot,
        }
    )


def save_results(
    output_path: str,
    track_result: BehaviorBootstrapResult,
    cell_result: BehaviorBootstrapResult,
    hierarchical_result: BehaviorBootstrapResult,
    save_estimates: bool = True,
) -> None:
    summary_track = result_to_dataframe(track_result)
    summary_cell = result_to_dataframe(cell_result)
    summary_hier = result_to_dataframe(hierarchical_result)

    summary_all = pd.concat(
        [summary_track, summary_cell, summary_hier],
        ignore_index=True,
    )

    summary_csv = os.path.join(output_path, "RESULT_bound_fraction_bootstrap.csv")
    summary_txt = os.path.join(output_path, "RESULT_bound_fraction_bootstrap.txt")

    summary_all.to_csv(summary_csv, index=False)

    with open(summary_txt, "w") as f:
        f.write("[Bound Fraction Bootstrapping]\n\n")
        f.write("Track-level bootstrap\n")
        f.write(summary_track.to_string(index=False))
        f.write("\n\nCell-level bootstrap\n")
        f.write(summary_cell.to_string(index=False))
        f.write("\n\nHierarchical bootstrap (cell -> track)\n")
        f.write(summary_hier.to_string(index=False))
        f.write("\n")

    if save_estimates:
        track_estimates_csv = os.path.join(output_path, "bound_fraction_bootstrap_track_estimates.csv")
        cell_estimates_csv = os.path.join(output_path, "bound_fraction_bootstrap_cell_estimates.csv")
        hier_estimates_csv = os.path.join(output_path, "bound_fraction_bootstrap_hierarchical_estimates.csv")

        track_result.bootstrap_estimates.to_csv(track_estimates_csv, index=False)
        cell_result.bootstrap_estimates.to_csv(cell_estimates_csv, index=False)
        hierarchical_result.bootstrap_estimates.to_csv(hier_estimates_csv, index=False)


def main(config_path: str = None) -> None:
    if not config_path:
        config_path = os.path.join(os.getcwd(), "script-config.toml")

    with open(config_path, "rb") as config_file:
        configs = tomllib.load(config_file)

    csv_path = configs["path"]["csv_path"]
    output_folder_name = configs["path"]["output_folder_name"]
    output_path = os.path.join(csv_path, output_folder_name)

    toggles = configs.get("toggle", {})
    use_gap_fixed = toggles.get("use_gap_fixed", False)

    bootstrap_cfg = configs.get("bound-fraction-bootstrap", {})
    n_boot = int(bootstrap_cfg.get("n_boot", 1000))
    random_state = bootstrap_cfg.get("random_state", 42)
    save_estimates = bool(bootstrap_cfg.get("save_bootstrap_estimates", True))

    logging_setup(output_path, "bound_fraction_bootstrap")

    input_csv = os.path.join(
        output_path,
        "intermediates",
        "gaps-and-fixes_decisions.csv" if use_gap_fixed else "bound_decisions.csv",
    )

    if not os.path.isfile(input_csv):
        raise FileNotFoundError(f"Input file not found: {input_csv}")

    print_log("[Bound Fraction Bootstrapping]")
    print_log("Reading:", input_csv)
    print_log("Bootstrap replicates:", n_boot)
    print_log("Random state:", random_state)
    print_log("Use gap-fixed input:", use_gap_fixed)

    spots_df = pd.read_csv(input_csv)
    spots_df = spots_df.loc[:, ~spots_df.columns.str.contains("^Unnamed")]

    track_result = bootstrap_behavior_proportions_by_track(
        spots_df=spots_df,
        n_boot=n_boot,
        random_state=random_state,
    )
    cell_result = bootstrap_behavior_proportions_by_cell(
        spots_df=spots_df,
        n_boot=n_boot,
        random_state=random_state,
    )
    hierarchical_result = bootstrap_behavior_proportions_hierarchical(
        spots_df=spots_df,
        n_boot=n_boot,
        random_state=random_state,
    )

    track_summary_df = result_to_dataframe(track_result)
    cell_summary_df = result_to_dataframe(cell_result)
    hier_summary_df = result_to_dataframe(hierarchical_result)

    print_log("\nTrack-level bootstrap")
    print_log(track_summary_df.to_string(index=False))

    print_log("\nCell-level bootstrap")
    print_log(cell_summary_df.to_string(index=False))

    print_log("\nHierarchical bootstrap (cell -> track)")
    print_log(hier_summary_df.to_string(index=False))

    save_results(
        output_path=output_path,
        track_result=track_result,
        cell_result=cell_result,
        hierarchical_result=hierarchical_result,
        save_estimates=save_estimates,
    )

    print_log("\nSaved:")
    print_log(os.path.join(output_path, "RESULT_bound_fraction_bootstrap.csv"))
    print_log(os.path.join(output_path, "RESULT_bound_fraction_bootstrap.txt"))
    if save_estimates:
        print_log(os.path.join(output_path, "bound_fraction_bootstrap_track_estimates.csv"))
        print_log(os.path.join(output_path, "bound_fraction_bootstrap_cell_estimates.csv"))
        print_log(os.path.join(output_path, "bound_fraction_bootstrap_hierarchical_estimates.csv"))


if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser(
        prog="bound_frac_bootstrapping",
        description="Bootstrap frame-level behavior proportions for SMol-FIESTA.",
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="Path to TOML config file.",
    )
    args = parser.parse_args()

    main(args.config)

    print(f"\n--- {time.time() - start_time:.2f} seconds ---")