# PyCBC Core Code Review: Multi-Detector Search & Asymmetric Windowing
**Branch:** `firinspiral3-multidet-asym`
**Commits:** `2dc84d3b` and `425c49c6`

## 1. Consistency with PyCBC Conventions
- **CLI Option Naming**: The addition of `--instruments` explicitly taking multiple space-separated strings (`nargs="+"`) conforms well with modern PyCBC and IGWN CLI patterns. The `--asym-*` options correctly use PyCBC naming paradigms (e.g., `--asym-windowing`, `--asym-snr-threshold`, etc.) inside the standard parser groups.
- **HDF5 Output Schema**: Writing single-detector features under their proper `f['H1/...']` and `f['L1/...']` hierarchies is well implemented via the `H5FileSyntSugar` wrapper, which intelligently switches to `mode='a'` allowing joint append operations to the same HDF5 file using an IFO prefix. The integration of `search/non_threshold_time` is also standard and follows the standard `search/*` metadata schema.
- **Post-processing Compatibility**: The `mergetrigs`, `findtrigs`, and `statmap` executables cleanly read and cascade `search/non_threshold_time`. Merging logic handles attribute aggregation safely, translating values appropriately to `all_trigs.attrs` for fair livetime estimation.

## 2. Minimality and Code Cleanliness
- **Modifications to Core Libraries**:
  - `pycbc/events/eventmgr.py`: Clean update to `H5FileSyntSugar` introducing a context manager and modifying behavior to support nested prefixes gracefully without sweeping re-architectures.
  - `pycbc/strain/strain.py`: `from_cli_multi_ifos` loops over `from_cli_single_ifo` using dictionary mappings. Minimal and functional.
  - `pycbc/workflow/matched_filter.py`: Properly implements a tiling and job submission generation routine (`setup_matchedfltr_dax_generated_multi`) supporting `--coincident-only` seglists.
- **Formatting and Cruft**: 
  - There are no latent `print()` statements or leftover debugging markers. 
  - Unused variable checking reveals tight variable scopes (e.g., `group_engine` dynamically mapped to `l1_engines` directly drops reference correctly). 

## 3. Algorithmic Correctness
- **`sample_asym_bins` (Asymmetric Windowing)**: 
  - Iterates effectively over 100 bins of 0.1s centered dynamically (`win_start_samp = t_center - half_window_samp`).
  - Block FFT overlaps seamlessly, caching via `block_f_cache` avoiding forward FFT duplication.
  - Properly extracts peak `sub_max_idx = np.argmax(np.abs(sub_corr))` independent of SNR thresholds, strictly following the un-thresholded local maximum requirement. Interval accumulation (`searched_intervals.append((gps_s, gps_e))`) appropriately logs GPS times for livetime aggregation downstream.
- **Fair Background Livetime**: 
  - Computed seamlessly in `statmap` as `(t_obs1 * t_nothresh2 + t_obs2 * t_nothresh1) / ts_interval`.
  - Maps to the theoretical formula $T_{\text{bkg}} = \frac{T_A \times T_{B,\text{no\_thresh}} + T_B \times T_{A,\text{no\_thresh}}}{\Delta t_{\text{slide}}}$. Algebraically robust mapping `t_obs1` to $T_A$ and `t_nothresh1` to $T_{A, \text{no\_thresh}}$. 

## 4. Memory and Performance
- **Joint Multi-detector Buffer Reuse**: Multi-IFO operations iterate (`for ifo in active_ifos:`) inside the hot segment loops in `pycbc_inspiral_fir` invoking `group_engine.process_segment(...)` while sharing the identical instantiated `group_engine` (from `MatchedFilterRatioControl`). Temporary FFT arrays (`corr_view`, `mult_view`) inside `group_engine._get_ifft_plan` are dimensionally cached by block size mapping, completely precluding double heap allocation for successive IFOs analyzing identical block structures.
- **`_ap_ref_key` Caching**: The engine caches using `(id(profile), id(psd))` in `MatchedFilterRatioControl`. Since ID pointers for references and static PSD arrays rarely shift within a single templated execution loop, it perfectly averts redundant recomputation and bank re-evaluations during hierarchical/coarse passes.

## Recommendation
The architecture implements the requested algorithmic requirements accurately. Code is optimally parallelized and aggressively circumvents cache misses and memory reallocation. Approve merge with no architectural or logical blockers.
