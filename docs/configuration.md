### Paths

| Parameter            | Type   | Description                                                        | Notes & use case                                                                 |
|----------------------|--------|--------------------------------------------------------------------|----------------------------------------------------------------------------------|
| `csv_path`           | string | Folder containing your TrackMate CSV exports.                      | Must contain `*_Cell_*_spots.csv`; if empty, scripts will error out immediately. |
| `mask_path`          | string | Folder of per-frame PNG masks; gray value = cell ID.               | Dimensions must match video frames—otherwise, overlays misalign.                 |
| `comdet_path`        | string | (Optional) Folder containing ComDet CSVs for fixed-particle analysis. | Leave blank if not using Fixed particle analysis                                 |
| `output_folder_name` | string | Name for the analysis output directory (created/overwritten).       | Changing it preserves old results but may break downstream lookups.              |

---

### Toggle Flags

| Parameter                 | Type | Description                                                                                  | Notes & use case                                                       |
|---------------------------|------|----------------------------------------------------------------------------------------------|------------------------------------------------------------------------|
| `use_gap_fixed`           | bool | Interpolate short detection gaps before reclassification.                                     | Compare runs with & without: flip this to gauge effect of gap-filling. |
| `filter_by_binding_prop`  | bool | Drop tracks whose track-level (bound+constrained) fraction is outside `[min_prop_binding, max_prop_binding]`. | If you see zero output, try disabling or widening these bounds.        |
| `run_visualizer`          | bool | Generate overlay videos showing track states (`green=bound`, `blue=constrained`, `red=fast`). | Disable for headless/batch runs to save time.                          |
| `run_MSD_calculations`    | bool | Compute & plot MSD curves for each track.                                                    | Turn off if you only need rebinding stats.                             |
| `run_fixed_particle`      | bool | Execute the (experimental) ComDet-based fixed-particle rebinding script.                     | Only effective when `comdet_path` is set.                              |
| `bacteria_analysis`       | bool | Aggressive overlap filtering, assuming one molecule per cell (ideal for bacterial data).     | Turn off for eukaryotic / high-density data.                           |

---

### Track-Sorting

| Parameter                   | Type | Description                                                                                       | Notes & use case                                                                                           |
|-----------------------------|------|---------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|
| `allowed_gap_max`           | int  | Max consecutive missing frames before splitting a track.                                          | Increase to tolerate long fluorophore dark states or noisy data; decrease to avoid fusing distinct events. |
| `allowed_track_length_min`  | int  | Discard tracks shorter than this (frames).                                                        | Raise to remove noisy, too-short tracks.                                                                   |
| `allowed_track_length_max`  | int  | Discard tracks longer than this (frames).                                                         | Lower to exclude immobile artefacts (e.g., stuck dirt).                                                    |
| `dist_range`                | int  | Frames before/after for computing inter-frame displacements.                                      | Larger windows → smoother classification but risk missing fast transitions.                                |
| `allowed_overlap`           | int  | Overlap-frames threshold for splitting/merging collisions (bacteria mode).                        |                                                                                                            |
| `concurrent_max`            | int  | Max consecutive frames with >1 spot before removing frames from the track (bacteria mode).       |                                                                                                            |

---

### Bound-Classification

| Parameter                   | Type   | Description                                                                      | Notes & use case                                              |
|-----------------------------|--------|----------------------------------------------------------------------------------|---------------------------------------------------------------|
| `distance_threshold`        | float  | px: relaxed step-size threshold over contextual frames to call “constrained” diffusion. | Most of the time, 2× `distance_threshold_strict` should work. |
| `distance_threshold_strict` | float  | px: strict step-size threshold over ±4-frame window to calling “bound”.          | Use the formula to find the best for your specific scenario.  |

---

### Gaps & Fixes

| Parameter               | Type   | Description                                                                    | Notes & use case                                                                             |
|-------------------------|--------|--------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
| `min_time_strict`       | int    | Min frames for strict-binding events; shorter ⇒ reclassify to constrained.     | Raise to suppress brief bindings.                                                            |
| `min_time_constrained`  | int    | Min frames for constrained diffusion; shorter ⇒ free diffusion.                 | Raise to suppress brief false constrained diffusion.                                         |
| `min_time_diffusion`    | int    | Min frames for free diffusion; shorter ⇒ context-based reclassification.       | Rarely > 0.                                                                                  |
| `max_bound_gapFill`     | int    | Max number of frames of gap allowed to be reclassified into binding events.    | Too high → merges distinct binds; too low → leaves genuine misses unfilled.                  |
| `min_prop_binding`      | float  | Min track-level bound+constrained fraction to retain a track.                  | Raise if you are testing populations of molecules that present binding or transient binding. |
| `max_prop_binding`      | float  | Max track-level bound+constrained fraction to retain a track.                  | Set > 0.99 to drop permanently stuck tracks or debris.                                       |

---

### Rebind-Analysis

| Parameter                         | Type   | Description                                                                                     | Notes & use case                                                                                    |
|-----------------------------------|--------|-------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------|
| `rebind_distance_same`            | float  | px: distance threshold for “same-site” rebinding (no fiducial).                                 | Set ≈ `distance_threshold_strict`.                                                                  |
| `rebind_distance_diff`            | float  | px: distance threshold for “different binding site” rebinding.                                  | Larger → a molecule needs to bind away for a rebinding to be counted as a different starting point. |
| `min_time_bound_strict`           | int    | Min frames per strict-binding to include in rebinding logs.                                     | Excludes fleeting contacts.                                                                         |
| `min_time_bound_constrained`      | int    | Min frames per constrained-binding to include in logs.                                           | Same rationale as strict, for the constrained class.                                                |                                       |
| `min_time_rebinding_strict`       | int    | Strict min diffusion frames to avoid counting flicker as rebinding.                             | Increase to filter out false positives.                                                             |
| `min_time_diffusion_subsequent`   | int    | Min diffusion frames after a bind to count as a “subsequent” event.                             | Ensures a real diffusion gap.                                                                       |
| `max_time_rebinding`              | int    | Frame window to search for rebinding; beyond this → unsuccessful.                                |                                                                                                     |

---

### MSD-Analysis

| Parameter          | Type   | Description                                                | Notes & use case                                         |
|--------------------|--------|------------------------------------------------------------|----------------------------------------------------------|
| `time_interval_ms` | float  | Frame interval in ms.                                       | Must exactly match the acquisition setting.              |
| `pixel_size_um`    | float  | Physical pixel size in µm.                                  | Critical for correct diffusion-coefficient calculations. |
| `use_brownian_only`| bool   | Force anomalous exponent α = 1 (pure Brownian fit).         | Disable to let the fit return both D and α.              |

---

### Visualizer

| Parameter             | Type   | Description                                                         | Notes & use case                                      |
|-----------------------|--------|---------------------------------------------------------------------|-------------------------------------------------------|
| `tail_effect_enabled` | bool   | Draw motion-blur “tails” for moving spots.                           | Disable if visual clutter or performance is an issue. |
| `spot_diameter`       | int    | Overlay spot diameter in pixels.                                     | Adjust to match the size of the real spot.            |

---

### Rebind-Fixed-Particle

| Parameter               | Type   | Description                                                                                   | Notes & use case                                                            |
|-------------------------|--------|-----------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------|
| `allowed_spot_overlap`  | float  | Max fractional area overlap for ComDet spots within a cell before skipping.                   | Lower to enforce strict spot isolation; raise if spot shapes vary slightly. |
