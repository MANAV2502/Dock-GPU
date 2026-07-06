# Dock-GPU: GPU-Accelerated Rigid/Flexible Receptor Docking Pipeline

**Notebook file:** `Dock_GPU.ipynb`
**Companion execution scripts:** `gpu_rigid_colab_docking.sh` · `gpu_flexible_colab_docking.sh`

---

## 1. Overview

Dock-GPU is a Google Colab–native front end for **AutoDock-GPU**, engineered to take a medicinal chemist or structural biologist from a prepared receptor/ligand set to ranked binding poses with zero local compute infrastructure and no command-line proficiency required. The notebook (`Dock_GPU.ipynb`) handles environment provisioning, file staging, and result retrieval through parameterized `#@param` cells, while the actual docking logic — GPU dispatch, per-ligand validation, pose extraction, and run reporting — is delegated to one of two production-grade Bash scripts, selected according to whether the target receptor is treated as **rigid** or carries **flexible side-chain residues**.

This separation of concerns is deliberate: the notebook is a *staging and orchestration layer*, while the shell scripts are the *docking engine*, independently versioned and independently auditable. A user can, in principle, run either script on any Linux machine with a CUDA-capable GPU outside of Colab entirely.

---

## 2. Repository Contents

| File | Role |
|---|---|
| `Dock_GPU.ipynb` | Colab notebook (Dock-GPU UI). Handles GPU/CUDA verification, AutoDock-GPU compilation, receptor/ligand/script upload, run orchestration, and results packaging. |
| `gpu_rigid_colab_docking.sh` | Docking engine for **rigid receptors** — no flexible side chains modeled during the search. |
| `gpu_flexible_colab_docking.sh` | Docking engine for **flexible-residue receptors** — auto-detects and injects `--flexres`, and additionally performs merged ligand+flexres best-pose extraction. |

Only one shell script is uploaded per session (Step 6 of the notebook); the notebook detects and executes whichever `.sh` file is present, so the same UI cells serve both docking modes.

---

## 3. Scientific Scope and Intended Use

This pipeline targets **structure-based virtual screening**: docking a library of prepared ligand `PDBQT` files against one or more prepared receptor grid maps (`AutoGrid4`-derived `.fld`/`.map` sets) using the GPU-parallelized Lamarckian Genetic Algorithm implemented in AutoDock-GPU. It assumes receptor and ligand preparation (protonation state assignment, partial charge calculation, torsion tree definition, grid box placement) has already been performed upstream — e.g., via ADFRsuite (`prepare_receptor4.py`, `prepare_flexreceptor4.py`) and a Meeko- or AutoDockTools-based ligand preparation workflow. Dock-GPU does not perform receptor/ligand preparation itself; it consumes their outputs.

For flexible-residue campaigns, the receptor is expected to be split into a **rigid core PDBQT** and a **flexible-residue PDBQT** (`<receptor_name>_flex.pdbqt`), consistent with the `prepare_flexreceptor4.py` convention — this naming pattern is what the flexible script's auto-detection logic keys on.

---

## 4. Prerequisites

Before opening the notebook, have the following ready, ideally pre-zipped for fast upload:

- **Receptor package(s):** AutoGrid `.fld` map file, all referenced `.map` files, and the rigid receptor `.pdbqt` (plus `<name>_flex.pdbqt` if docking flexibly). Each receptor should occupy its own named subfolder if multiple targets are screened in one session.
- **Ligand library:** individual `.pdbqt` files, one per compound, already protonated and torsion-tree-assigned.
- **Docking script:** either `gpu_rigid_colab_docking.sh` or `gpu_flexible_colab_docking.sh`, matching the receptor preparation mode.
- **A Colab GPU runtime** (Runtime → Change runtime type → GPU). The notebook's Step 1 cells verify GPU presence, CUDA toolkit version, and GCC availability — the compiled CUDA version **must** match the runtime's `nvcc` version, or the AutoDock-GPU build in Step 2 will fail or silently mis-target the wrong compute capability.

---

## 5. Notebook Walkthrough (`Dock_GPU.ipynb`)

The notebook is organized as ten sequential, self-documenting steps. Each executable cell is annotated with a `#@title` and `#@markdown` description visible directly in the Colab UI.

| Step | Purpose |
|---|---|
| **1 — System & GPU Verification** | Confirms GPU allocation (`nvidia-smi`), CUDA compiler version (`nvcc --version`), and GCC availability. The noted CUDA version must be carried into Step 2. |
| **2 — Install & Compile AutoDock-GPU** | Clones `ccsb-scripps/AutoDock-GPU` from source and compiles with `DEVICE=CUDA`, `NUMWI=256`, matching the CUDA path detected in Step 1. |
| **3 — Directory Setup** | Creates and enters the working `docking/` directory that anchors all subsequent relative paths. |
| **4 — Receptor Upload** | Prompts for a `protein_name` and one of three upload modes (individual file picker, ZIP archive, or Google Drive). Deposits files flat into `docking/rigid_flexible/<protein_name>/`, supporting multiple named receptors within a single session without collision. |
| **5 — Ligand Upload** | Same three upload modes, depositing `.pdbqt` ligands flat into `docking/ligand/`. |
| **6 — Docking Script Upload** | Accepts either shell script under any filename; the notebook persists the chosen filename to `.script_name` so later cells recover it correctly even after a runtime restart. |
| **7 — File-Count Verification** | Cross-checks ligand and per-receptor file counts before committing GPU time to a run. |
| **8 — Run AutoDock-GPU** | Normalizes line endings (CRLF → LF), sets execute permissions, and launches the uploaded script. This is the long-running compute step. |
| **9 — Verify & Download Results** | Confirms output file counts (expected: 3 files per ligand — `.dlg`, `.xml`, best-pose `.pdbqt`) per receptor, then archives everything to `results.tar.gz` for local download. |
| **10 — Batch Reset (optional)** | Clears the ligand directory and prior results to stage a second ligand batch against the same receptor without re-uploading receptor/map files. |

**Design notes worth flagging to a collaborator inheriting this notebook:**
- The receptor- and ligand-upload cells use a `_move_flat()` / `_extract_flat()` pattern that strips directory nesting on extraction — this is intentional to guarantee AutoDock-GPU's flat-directory discovery assumptions hold regardless of how a user's ZIP was structured locally.
- Script-name persistence (`.script_name` file) is what allows Steps 8 onward to survive a Colab runtime disconnect without forcing a re-upload of the docking script.

---

## 6. Docking Engine: Shared Architecture

Both `gpu_rigid_colab_docking.sh` and `gpu_flexible_colab_docking.sh` share a common nine-section chassis:

1. **User-configurable paths** — `BIN_DIR`, `RECEPTOR_ROOT_DIR`, `LIGAND_DIR`, `OUTPUT_ROOT`.
2. **Docking parameters** — passed directly to the `autodock_gpu_<N>wi` binary.
3. **Script behaviour controls** — timeout, energy sanity threshold, minimum free disk space.
4. **Output/log initialization** — timestamped log and CSV summary, with `tee` mirroring to stdout.
5. **Helper functions** — FLD map-integrity validation, ligand PDBQT keyword validation, best-energy extraction from DLG output, and an energy sanity filter for physically implausible binding energies.
6. **Pre-flight checks** — GPU/CUDA availability, auto-selection of the correct `autodock_gpu_{64,128,256}wi` binary, version logging for reproducibility, directory existence, and disk space.
7. **Ligand pre-validation** — every ligand PDBQT is validated once, up front, so malformed files are excluded before entering the docking loop rather than failing mid-run.
8. **Docking loop** — iterates every receptor × ligand pair, dispatches to the AutoDock-GPU binary, and records per-job status to the CSV.
9. **Final summary** — aggregate counts of successes, failures, timeouts, and (for the flexible script) flexres-enabled receptors and poses extracted.

### 6.1 Default Docking Parameters (Section 2)

| Flag | Default | Meaning |
|---|---|---|
| `--nev` | 250000 | Maximum number of energy evaluations per run |
| `--heurmax` | 1000000 | Ceiling for the heuristics-based evaluation estimate |
| `--nrun` | 5 | Independent LGA runs per ligand |
| `-p` (population) | 75 | GA population size |
| `-g` (generations) | 5000 | Maximum generations per run |
| `--lsmet` | `ad` (ADADELTA) | Local search method |
| `--lsit` | 750 | Local search iterations |
| `--autostop` | 1 (enabled) | Early termination on convergence |
| `--asfreq` | 5 | Autostop check frequency |
| `--stopstd` | 0.1 | Convergence standard-deviation threshold |
| `--clustering` | 1 (enabled) | Cluster docked poses by RMSD |
| `--gbest` | 1 (enabled) | Track global best pose across runs |
| `--xmloutput` / `--dlgoutput` | 1 / 1 | Emit XML and DLG result files |

Both `TIMEOUT_SECONDS` (default 300 s per ligand) and `ENERGY_SANITY_THRESHOLD` (default −20.0 kcal/mol, flagging suspiciously favorable predicted affinities that often indicate a grid or protonation artifact rather than genuine high-affinity binding) are exposed for tuning per campaign.

### 6.2 What Distinguishes the Flexible Script

`gpu_flexible_colab_docking.sh` extends the shared chassis with a fourth integration tier, tagged `[R4]` in-line:

- **Universal `--flexres` auto-detection.** For each receptor subfolder, the script looks for `<receptor_name>_flex.pdbqt`. If found and structurally valid, it is injected via `--flexres` automatically — no manual per-receptor flag editing required, even across a multi-receptor batch.
- **Flexible-residue PDBQT structural validation.** Unlike ligand PDBQTs, flexres files carry no global `TORSDOF`; each residue instead nests its own `ROOT`/`BRANCH`/`ENDBRANCH` inside `BEGIN_RES`/`END_RES` blocks. `validate_flexres_pdbqt()` checks this structure explicitly before the file is trusted.
- **Merged best-pose extraction.** Because AutoDock-GPU writes ligand atoms and co-evaluated flexible-residue atoms into the *same* `DOCKED: MODEL...ENDMDL` block, the script extracts the single lowest-energy run (identical selection criterion to the CSV-reported best energy) and writes it as a standalone PDBQT containing both ligand and flexres atoms — ready for direct MD system building without a separate stitching step.
- **Atom-count QC cross-check.** The extracted pose's ligand and flexres atom counts are cross-validated against the "Number of flexres atoms" line reported in the DLG header, flagging any silent truncation.
- **Graceful `python3` degradation.** Best-pose extraction relies on `python3` for DLG parsing; if unavailable, docking proceeds normally and only that one post-processing step is skipped (logged as a warning, not a failure).
- **Extended CSV schema.** The flexible script's summary CSV adds `flexres_used` and `best_pose_pdbqt` columns alongside the shared `receptor, ligand, status, best_energy_kcal_mol, dlg_path, notes` fields.

---

## 7. Output Structure

After a completed run, `docking/results_autodock_gpu/<receptor_name>/` contains, per ligand:

- `<ligand>.dlg` — full docking log (all runs, clustering, binding energies)
- `<ligand>.xml` — machine-readable run summary
- `<ligand>_best.pdbqt` — extracted lowest-energy pose (ligand-only for rigid runs; ligand+flexres merged for flexible runs)

A run-level `docking_summary_<timestamp>.csv` and `docking_run_<timestamp>.log` are written to `results_autodock_gpu/` for downstream triage and reproducibility record-keeping. All of this is bundled into `results.tar.gz` by the notebook's Step 9 for local download.

---

## 8. Troubleshooting Notes

- **Compilation failure in Step 2:** almost always a CUDA version mismatch between the Colab runtime's `nvcc` and the hardcoded `cuda-12.8` paths in the compile cell — update both occurrences to match Step 1's reported version.
- **`skipped_no_fld` / `skipped_ambiguous_fld` in the CSV:** the receptor subfolder is missing its `.maps.fld` file, or contains more than one — each receptor subfolder must resolve to exactly one FLD.
- **`--flexres` silently not applied:** confirm the flexible-residue file is named exactly `<receptor_name>_flex.pdbqt` and sits alongside the rigid PDBQT in the same subfolder; the auto-detection is filename-pattern-based, not content-based.
- **Result count ≠ 3× ligand count:** check the per-ligand status column in the summary CSV for `timeout`, `failed`, `no_dlg`, or `no_poses` entries before assuming a systemic pipeline fault.

---

## 9. Attribution

Docking is performed by **AutoDock-GPU**, sourced from `ccsb-scripps/AutoDock-GPU`. This repository provides only the orchestration notebook and production-hardened execution scripts around that engine; the underlying scoring function, search algorithm, and grid-map format are unmodified from upstream AutoDock-GPU. Receptor flexibility handling (`--flexres`) rests on the Flexibility Tree formalism introduced by AutoDockFR and implemented in ADFRsuite's `prepare_flexreceptor4.py`.

---

## 10. Citations

If results generated with this pipeline are used in a publication, thesis, or grant report, the underlying methods should be cited as follows.

**AutoDock-GPU (primary docking engine):**
Santos-Martins, D., Solis-Vasquez, L., Tillack, A. F., Sanner, M. F., Koch, A., & Forli, S. (2021). Accelerating AutoDock4 with GPUs and gradient-based local search. *Journal of Chemical Theory and Computation*, 17(2), 1060–1073. https://doi.org/10.1021/acs.jctc.0c01006

**AutoDock4 scoring function (foundation of the AutoDock-GPU force field):**
Morris, G. M., Huey, R., Lindstrom, W., Sanner, M. F., Belew, R. K., Goodsell, D. S., & Olson, A. J. (2009). AutoDock4 and AutoDockTools4: Automated docking with selective receptor flexibility. *Journal of Computational Chemistry*, 30(16), 2785–2791. https://doi.org/10.1002/jcc.21256

**AutoDockFR / ADFRsuite (flexible-residue receptor preparation — `prepare_receptor4.py`, `prepare_flexreceptor4.py`):**
Ravindranath, P. A., Forli, S., Goodsell, D. S., Olson, A. J., & Sanner, M. F. (2015). AutoDockFR: Advances in protein-ligand docking with explicitly specified binding site flexibility. *PLOS Computational Biology*, 11(12), e1004586. https://doi.org/10.1371/journal.pcbi.1004586

**Meeko (ligand/receptor PDBQT preparation toolkit, Forli lab):**
Meeko: Python package for preparing small molecules and receptors for AutoDock. Forli Lab, Scripps Research. Source and documentation: https://github.com/forlilab/Meeko

**GPU-accelerated large-scale deployment context (optional, if benchmarking or scale-up methodology is referenced):**
LeGrand, S., Scheinberg, A., Tillack, A. F., Thavappiragasam, M., Vermaas, J. V., Agarwal, R., Larkin, J., Poole, D., Santos-Martins, D., Solis-Vasquez, L., Koch, A., Forli, S., Hernandez, O., Smith, J. C., & Sedova, A. (2020). GPU-accelerated drug discovery with docking on the Summit supercomputer: Porting, optimization, and application to COVID-19 research. In *Proceedings of the 11th ACM International Conference on Bioinformatics, Computational Biology and Health Informatics* (BCB '20).

**Software repositories (for reproducibility/version pinning):**
- AutoDock-GPU: https://github.com/ccsb-scripps/AutoDock-GPU
- ADFRsuite: https://ccsb.scripps.edu/adfr/
- Meeko: https://github.com/forlilab/Meeko

When reporting results, cite the specific AutoDock-GPU release/commit hash compiled in Step 2 of the notebook, since scoring and performance characteristics can shift meaningfully across versions.
