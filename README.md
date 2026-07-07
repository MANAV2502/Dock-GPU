# Dock-GPU: GPU-Accelerated Rigid/Flexible Receptor Docking Pipeline

**Notebook:** `Dock_GPU.ipynb`

**Engines:** `gpu_rigid_colab_docking.sh` · `gpu_flexible_colab_docking.sh`

## Overview

Dock-GPU is a Colab front end for **AutoDock-GPU**. The notebook handles GPU/CUDA setup, file staging, and results retrieval via `#@param` cells; one of two Bash scripts (`gpu_rigid_colab_docking.sh`/`gpu_flexible_colab_docking.sh`) performs the actual docking, chosen by whether the receptor is rigid or carries flexible side chains. Either script can also run standalone on any CUDA-capable Linux machine.

**Scope:** Structure-based virtual screening, docking prepared ligand `.pdbqt` files against `AutoGrid4` receptor maps via AutoDock-GPU's GPU-parallelized LGA. Receptor/ligand preparation (protonation, charges, torsion trees, grid boxes) must already be done upstream. For flexible campaigns, the receptor must be split into a rigid PDBQT and a `<receptor_name>_flex.pdbqt`, this is what the flexible script's auto-detection keys on. Can be used for more than one receptor at a time, i.e, docking same set of ligand for more than one receptors.

## Prerequisites

- **Receptor:** `.fld` + `.map` files and rigid `.pdbqt` (plus `_flex.pdbqt` if flexible), one named subfolder per receptor
- **Ligands:** individual `.pdbqt` files, already protonated with torsion trees assigned
- **Docking script:** matching rigid/flexible mode (`gpu_rigid_colab_docking.sh`/`gpu_flexible_colab_docking.sh`)
- **Colab GPU runtime:** compiled AutoDock-GPU's CUDA version must match the runtime's `nvcc` version (checked in Step 1/2)

## Notebook Steps

| Step | Purpose |
|---|---|
| 1. System check | GPU, CUDA, GCC verification |
| 2. Build | Clone + compile AutoDock-GPU (`DEVICE=CUDA`, `NUMWI=256`) |
| 3. Directory setup | Creates `docking/` working directory |
| 4. Receptor upload | File picker / ZIP / Drive → `docking/rigid_flexible/<protein_name>/` |
| 5. Ligand upload | Same options → `docking/ligand/` |
| 6. Script upload | Any `.sh` filename accepted; name persisted for later steps |
| 7. File-count check | Confirms uploads before committing GPU time |
| 8. Run docking | Fixes line endings, sets permissions, launches script |
| 9. Results | Verifies output counts (3 files/ligand), packages `results.tar.gz` |
|10. Batch reset | Clears ligand/results dirs to re-run against the same receptor |

Uploads are flattened on extraction (no nested subfolders) so file discovery stays predictable regardless of local ZIP structure; the persisted script name lets Steps 8+ survive a runtime restart.

## Docking Script Architecture

Both scripts share: configurable paths → docking parameters → behaviour controls (timeout, energy sanity threshold, disk check) → logging/CSV init → validation helpers (FLD integrity, ligand PDBQT checks, best-energy extraction) → pre-flight checks (GPU, binary auto-select, versioning) → ligand pre-validation → docking loop → summary.

**Default parameters:** 
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

**Flexible script adds:**
- Auto-detects and injects `--flexres` from `<receptor_name>_flex.pdbqt`
- Validates flexres PDBQT structure (`BEGIN_RES`/`END_RES` blocks, no global `TORSDOF`)
- Extracts the merged ligand+flexres best pose directly from the shared `DOCKED:` block (no separate stitching step)
- Cross-checks extracted atom counts against the DLG header; degrades gracefully if `python3` is unavailable
- Extended CSV: adds `flexres_used`, `best_pose_pdbqt` columns

## Output

Per ligand, per receptor: `.dlg` (full log), `.xml` (machine-readable), `_best.pdbqt` (lowest-energy pose, ligand+flexres merged where applicable). A summary CSV and run log accompany each session, all bundled into `results.tar.gz`.

## Troubleshooting

- **Build fails:** CUDA version in the compile cell doesn't match `nvcc --version` — update both such that both are same
- **`skipped_no_fld`/`skipped_ambiguous_fld`:** receptor subfolder has zero or multiple `.maps.fld` files
- **`--flexres` not applied:** flexres file must be named exactly `<receptor_name>_flex.pdbqt` (detection is filename-based, not content-based)
- **Result count ≠ 3× ligands:** check the CSV status column for `timeout`/`failed`/`no_dlg`/`no_poses`

## Citations

- Santos-Martins, D. et al. (2021). Accelerating AutoDock4 with GPUs and gradient-based local search. *J. Chem. Theory Comput.*, 17(2), 1060–1073. https://doi.org/10.1021/acs.jctc.0c01006
- Morris, G. M. et al. (2009). AutoDock4 and AutoDockTools4: Automated docking with selective receptor flexibility. *J. Comput. Chem.*, 30(16), 2785–2791. https://doi.org/10.1002/jcc.21256
- Ravindranath, P. A. et al. (2015). AutoDockFR: Advances in protein-ligand docking with explicitly specified binding site flexibility. *PLOS Comput. Biol.*, 11(12), e1004586. https://doi.org/10.1371/journal.pcbi.1004586
- Meeko (Forli Lab, Scripps Research): https://github.com/forlilab/Meeko
- Repositories: AutoDock-GPU — https://github.com/ccsb-scripps/AutoDock-GPU · ADFRsuite — https://ccsb.scripps.edu/adfr/

Cite the specific AutoDock-GPU commit/release compiled in Step 2 — scoring behaviour is not guaranteed stable across versions.

## Citing This Pipeline

If Dock-GPU is used to generate docking results reported in a publication, thesis, or preprint, cite the repository directly, in addition to the underlying AutoDock-GPU/ADFRsuite/Meeko citations listed in `README.md`:

**Repository citation:**
Dock-GPU: A GPU-Accelerated Rigid/Flexible Receptor Docking Pipeline for Google Colab. GitHub repository: https://github.com/MANAV2502/Dock-GPU

**BibTeX:**
```bibtex
@software{dockgpu2026,
  title  = {Dock-GPU: A GPU-Accelerated Rigid/Flexible Receptor Docking Pipeline for Google Colab},
  author = {{Manav Patel}},
  year   = {2026},
  url    = {https://github.com/MANAV2502/Dock-GPU}
}
```

Include the specific commit hash of the repository used at the time of the docking run (e.g., via `git rev-parse HEAD` or the commit visible in the GitHub UI), since the rigid/flexible scripts and notebook may be revised over time and downstream reproducibility depends on pinning the exact version consulted.
