# chik_mapping

R project estimating global chikungunya (CHIKV) burden: force-of-infection (FOI) mapping, age-stratified infection/burden estimation, probabilistic sensitivity analysis (PSA), and relative-risk (RR) adjustment for multimorbidity, feeding into burden maps and at-risk population estimates.

## Layout

- Root: standalone analysis scripts (`chikmap_psa_final.R`, `foi_map.R`, `age_strat_burden*.R`, `cookie_cut_map.R`, `SpatialBlockBootstrap_v3.R`, RR-integration scripts, etc.) plus explanatory markdown docs (`ANALYSIS_*.md`, `IMPLEMENTATION_*.md`, `RR_INTEGRATION_GUIDE.md`).
- `Functions/` — shared function libraries (`BurdenFunctions.R` / `BurdenFunctions_v2.R`, thinning, raster plotting, NA-fixing helpers) sourced by the analysis scripts.
- `Scripts/` — a second, overlapping set of pipeline scripts (some names duplicate root-level scripts, e.g. `age_strat_burden_estim.R`, `cookie_cut_map.R`, `covariates_all.R`, `lhs_samples.R` — these are separate copies, not the same file; check dates/content before assuming which is current).
- `Data/` — small tracked CSVs (`global_all_atrisk.csv`, `global_all_focal.csv`).
- `MainData/` — large model inputs/outputs (FOI samples, combined burden estimates), several 1–3 GB `.RData` files. **Not tracked in git** (see below).

## Git conventions

- `*.RData` / `*.rds` anywhere in the repo are gitignored — they're multi-GB model outputs that exceed GitHub's file size limits. Share these via OneDrive/external storage, not git.
- Origin: `https://github.com/hyolimkang/chik_mapping`, default branch `main`.
- `gh` CLI is installed locally; run `gh auth login` once to enable PR/issue workflows from the terminal.

## Working conventions

- Run `/code-review` (or `/code-review ultra` for a deeper multi-agent pass) before committing non-trivial analysis changes.
- Run `/security-review` if a change touches credentials, external data fetches, or file/path handling.
- For visual outputs (burden maps, PSA summary dashboards), prefer building an Artifact (HTML) over a static plot when an interactive/shareable view is useful.
