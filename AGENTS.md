# AGENTS.md

Guidance for AI agents working in this project.

## Scope

This guide is scoped to the workflow defined by `exe.sh`. Treat `exe.sh` as the source of truth for how the project is built, fitted, and plotted. Do not duplicate its command sequence in notes or answers; read the current file when you need exact commands.

## What The Code Does

This is a ROOT/C++ and Python project for CMS jet calibration studies. The core executable reads a JSON run card, opens ROOT input files containing response and PF-composition graphs, applies FSR and systematic variations, runs a global JES/PF-composition fit with Minuit, writes fit products to a ROOT output file, and then a Python plotting script makes plots from that output.

At a high level, the project has three moving parts:

1. Build and run the C++ `GlobalFit` executable.
2. Configure the fit through a JSON run card and C++ maps.
3. Plot the ROOT output with the plotting script selected by `exe.sh`.

## Where To Look First

Start with these files, in this order:

1. `exe.sh`: authoritative build, fit, and plotting workflow.
2. `jsons/info.json`: run card passed to `GlobalFit`.
3. `Makefile`: executable build rules, compiler flags, ROOT flags, and linked libraries.
4. `src/GlobalFit_exe.cc`: executable wrapper that constructs `GlobalFit` and calls `Run()`.
5. `include/GlobalFit.hpp`: public fit pipeline, configuration constants, and global state.
6. `src/GlobalFit.cxx`: main implementation, data loading, fit setup, Minuit objective, and output writing.
7. `include/constants.hpp`: ROOT object maps, input names, nuisance maps, FSR maps, fit ranges, and analytic shape formulas.
8. `include/Containers.hpp` and `src/Containers.cxx`: data, shape, and nuisance container classes.
9. Plotting script named by `exe.sh`: ROOT output expectations and plotting choices.

## Fit Configuration

The JSON file selected by `exe.sh` controls the active fit. The current project layout uses `jsons/info.json`, but always check `exe.sh` before assuming that.

Important run-card fields:

- `run`, `mode`, `eta_min`, `eta_max`: replace `RUN`, `MODE`, `ETAMIN`, and `ETAMAX` placeholders in strings from `include/constants.hpp`.
- `samples`, `hdm_methods`, `types`: select which input graph names are activated.
- `shapes`: selects shape families from `shapes_map`.
- `output_fname`: names the ROOT output file. Keep the plotting script aligned with this value.

Input ROOT filename patterns are defined in `include/constants.hpp`; check `input_fnames` before assuming which input file is read.

## GlobalFit Pipeline

The executable path enters through `src/GlobalFit_exe.cc`, then runs `GlobalFit::Run()` in `src/GlobalFit.cxx`.

The pipeline is:

1. Load selected input graphs from `input_hnames_map`.
2. Load nuisance histograms from `nuisances_map`.
3. Build analytic shape functions from `shapes_map`.
4. Load reference histograms.
5. Load FSR correction histograms from `kfsr_hnames_map`.
6. Apply FSR shifts, including multijet recoil handling.
7. Apply graph range selections.
8. Build the global fit function and `TFitter`.
9. Run Minuit and retrieve fit summaries plus the error matrix.
10. Write fit objects, data graphs, variations, nuisance shapes, references, and the error matrix.

The Minuit objective is `jesFitter()` in `src/GlobalFit.cxx`. It loops over data points, applies nuisance shifts, compares against `_jesFit`, and adds Gaussian penalties for nuisance parameters and, when enabled, fit parameters.

## Data Model

`DataContainer` owns cloned graphs:

- `raw`: original graph after zero-point cleanup.
- `input`: prefit graph after corrections.
- `output`: postfit or nuisance-shifted graph.
- `variation`: graph used for uncertainty propagation.

`ShapeContainer` owns a `TF1` analytic shape and a fit-parameter index. Shape families can share one parameter across observables.

`NuisanceContainer` owns a cloned `TH1D` systematic variation and records which input graph it applies to.

## Plotting

Use the plotting script named by `exe.sh`. Check the script itself for:

- Which ROOT file it opens.
- Which output object names it expects from `StoreFitOutput()`.
- Which plot suffixes and style choices are set in the script.
- Whether missing ROOT objects are validated before drawing.

Keep the plotting script aligned with the run card's `output_fname` and the objects written by `StoreFitOutput()`.

## Common Edits

To add or change a global-fit input:

1. Add or adjust the ROOT object mapping in `include/constants.hpp`.
2. Ensure the name can be selected by the run-card fields.
3. Verify the object exists in the input ROOT file selected by `input_fnames` and the run card.
4. Update the plotting script if the output object should be plotted.

To add a new shape family:

1. Add base and observable-specific entries to `shapes_map` in `include/constants.hpp`.
2. Add the family name to the run card selected by `exe.sh`.
3. Add plotting entries if the shape should appear in plots.

To change run, eta bin, mode, input file, or output file:

1. Prefer editing the run card and central maps.
2. Keep the plotting script aligned with output file and object names.
3. Re-read `exe.sh` before running anything.

## Conventions And Sharp Edges

- `GlobalFit::current_obs` controls how `_jesFit` evaluates. Set it before evaluating the fit function manually.
- Response fits may use the reference JES histogram through `useJESref`; check the setting before interpreting output.
- PF-composition observables are scaled during fitting through `ScaleFullSimShape`; check the constant before changing plots or units.
- `jesFitter()` mutates `output` and `variation` graphs while evaluating the fit. Treat those graphs as final only after the global fit completes.
- The code uses global maps declared in `include/GlobalFit.hpp`: `shapes`, `my_data`, `nuisances`, and `recoils`. Do not assume multiple independent fits can run in one process without cleanup.
- Missing ROOT objects often fail by `assert`; builds with assertions disabled may hide broken inputs.
- `include/constants.hpp` is large. Search by exact shape or input name before editing.
- The plotting script may receive missing ROOT objects as `None`; check validation before drawing.
- Review `src/utils.cxx::FuncToHist()` before relying on histogram-valued fit output.

## Search Anchors

Useful symbols and maps:

- `GlobalFit::Run`
- `jesFitter`
- `jesFit_wrapper`
- `input_fnames`
- `input_hnames_map`
- `nuisances_map`
- `kfsr_hnames_map`
- `shapes_map`
- `output_fname`
- `Resp_comb`
