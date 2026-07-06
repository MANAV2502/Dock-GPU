#!/bin/bash
set -uo pipefail

############################################################
#  AutoDock-GPU Rigid Docking — Production Script
#  Google Colab | Universal + Precise + Resilient
#
#  Integrations (in order of addition):
#  [R1] Universal .maps.fld discovery (any filename)
#  [R2] Auto-binary selection, GPU check, ligand/receptor
#       discovery, disk check, counters, timestamped log
#  [R3] FLD map-integrity check, PDBQT validation, DLG
#       output verification, resume support, energy sanity
#       filter, per-ligand timeout, version logging,
#       CSV run report, receptor PDBQT presence warning
#  [R4] Universal --flexres auto-detection (any receptor_
#       name_flex.pdbqt), structural validation, and best-
#       pose extraction: pulls the single lowest-energy Run
#       out of the DLG and writes it as a standalone PDBQT
#       containing the ligand atoms AND the co-evaluated
#       flexible-residue atoms for that exact pose, with
#       atom-count QC cross-checks against the DLG header.
############################################################


# ═══════════════════════════════════════════════════════════
# SECTION 1 — USER-CONFIGURABLE PATHS
#   Only edit values in this section.
# ═══════════════════════════════════════════════════════════
BIN_DIR="/content/AutoDock-GPU/bin"
RECEPTOR_ROOT_DIR="/content/docking/rigid_flexible"
LIGAND_DIR="/content/docking/ligand"
OUTPUT_ROOT="/content/docking/results_autodock_gpu"


# ═══════════════════════════════════════════════════════════
# SECTION 2 — DOCKING PARAMETERS
# ═══════════════════════════════════════════════════════════
NEV="250000"
HEURMAX="1000000"
NRUN="5"
P="75"
G="5000"
LSMET="ad"
LSIT="750"
AUTOSTOP="1"
ASFREQ="5"
STOPSTD="0.1"
CLUSTERING="1"
GBEST="1"
XML="1"
DLG="1"


# ═══════════════════════════════════════════════════════════
# SECTION 3 — SCRIPT BEHAVIOUR CONTROLS
# ═══════════════════════════════════════════════════════════
TIMEOUT_SECONDS=300            # max seconds per ligand docking job
ENERGY_SANITY_THRESHOLD="-20.0"  # flag best energies below this (kcal/mol)
MIN_DISK_GB=1                 # minimum free disk space required (GB)


# ═══════════════════════════════════════════════════════════
# SECTION 4 — OUTPUT / LOG INITIALISATION
# ═══════════════════════════════════════════════════════════
mkdir -p "$OUTPUT_ROOT"
LOG="$OUTPUT_ROOT/docking_run_$(date +%Y%m%d_%H%M%S).log"
REPORT_CSV="$OUTPUT_ROOT/docking_summary_$(date +%Y%m%d_%H%M%S).csv"

# Mirror everything (stdout + stderr) to the log file
exec > >(tee -a "$LOG") 2>&1

echo "══════════════════════════════════════════════════════"
echo "  AutoDock-GPU Rigid_Flexible Docking — Production Script"
echo "  Started : $(date)"
echo "══════════════════════════════════════════════════════"
echo "  Receptor root : $RECEPTOR_ROOT_DIR"
echo "  Ligand dir    : $LIGAND_DIR"
echo "  Output root   : $OUTPUT_ROOT"
echo "  Log file      : $LOG"
echo "  CSV report    : $REPORT_CSV"
echo

# Initialise CSV with header
echo "receptor,ligand,status,best_energy_kcal_mol,dlg_path,notes,flexres_used,best_pose_pdbqt" > "$REPORT_CSV"


# ═══════════════════════════════════════════════════════════
# SECTION 5 — HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════

# [R3] Cross-validate every .map file referenced inside the FLD
validate_fld() {
    local fld="$1"
    local fld_dir
    fld_dir="$(dirname "$fld")"
    local missing=0

    while IFS= read -r line; do
        # FLD lines referencing maps look like: "map <receptor>.X.map"
        if [[ "$line" =~ ^map[[:space:]]+(.+\.map)$ ]]; then
            local map_file="$fld_dir/${BASH_REMATCH[1]}"
            if [[ ! -f "$map_file" ]]; then
                echo "     ✖ Missing map file: ${BASH_REMATCH[1]}"
                missing=$(( missing + 1 ))
            fi
        fi
    done < "$fld"

    [[ $missing -eq 0 ]]
}

# [R3] Check that a ligand PDBQT contains required AutoDock4 keywords
validate_pdbqt() {
    local pdbqt="$1"
    local errors=0

    grep -q "^ROOT"            "$pdbqt" || { echo "     ✖ Missing ROOT keyword";     errors=$(( errors + 1 )); }
    grep -q "^ENDROOT"         "$pdbqt" || { echo "     ✖ Missing ENDROOT keyword";  errors=$(( errors + 1 )); }
    grep -q "^TORSDOF"         "$pdbqt" || { echo "     ✖ Missing TORSDOF keyword";  errors=$(( errors + 1 )); }
    grep -qE "^ATOM|^HETATM"  "$pdbqt" || { echo "     ✖ No ATOM/HETATM records";   errors=$(( errors + 1 )); }

    [[ $errors -eq 0 ]]
}

# [R3] Extract the best (lowest) binding energy from a DLG file
get_best_energy() {
    local dlg="$1"
    grep "Free Energy of Binding" "$dlg" \
        | grep -oP '[+-]?[0-9]+\.[0-9]+' \
        | sort -n | head -1 2>/dev/null \
    || echo "N/A"
}

# [R3] Flag physically unrealistic binding energies
check_energy_sanity() {
    local energy="$1"
    local ligand_name="$2"
    if [[ "$energy" != "N/A" ]]; then
        if awk "BEGIN { exit !($energy < $ENERGY_SANITY_THRESHOLD) }"; then
            echo "   ⚠  Sanity warning : Best energy = ${energy} kcal/mol"
            echo "      Values below ${ENERGY_SANITY_THRESHOLD} kcal/mol may indicate"
            echo "      grid box misconfiguration or atom-type assignment errors."
        fi
    fi
}

# [R4] Validate flexible-residue PDBQT structure (prepare_flexreceptor4.py output).
#      Unlike ligand PDBQTs, flexres files have no global TORSDOF — each residue
#      carries its own ROOT/BRANCH/ENDBRANCH nest wrapped in BEGIN_RES/END_RES.
validate_flexres_pdbqt() {
    local pdbqt="$1"
    local errors=0

    grep -q "^BEGIN_RES"      "$pdbqt" || { echo "     ✖ Missing BEGIN_RES keyword"; errors=$(( errors + 1 )); }
    grep -q "^END_RES"        "$pdbqt" || { echo "     ✖ Missing END_RES keyword";   errors=$(( errors + 1 )); }
    grep -q "^ROOT"            "$pdbqt" || { echo "     ✖ Missing ROOT keyword";      errors=$(( errors + 1 )); }
    grep -qE "^ATOM|^HETATM"  "$pdbqt" || { echo "     ✖ No ATOM/HETATM records";    errors=$(( errors + 1 )); }

    # Unbalanced BEGIN_RES/END_RES markers indicate truncated or corrupted output
    local n_begin n_end
    n_begin=$(grep -c "^BEGIN_RES" "$pdbqt")
    n_end=$(grep -c "^END_RES"   "$pdbqt")
    if [[ "$n_begin" -ne "$n_end" ]]; then
        echo "     ✖ Unbalanced BEGIN_RES ($n_begin) / END_RES ($n_end) markers"
        errors=$(( errors + 1 ))
    fi

    [[ $errors -eq 0 ]]
}

# [R4] Extract the single lowest-energy pose from a DLG into a standalone PDBQT.
#      Uses the identical "lowest Estimated Free Energy of Binding across all
#      Runs" criterion as get_best_energy() above, so the extracted pose always
#      corresponds to the same Run already reported in the CSV. When --flexres
#      was used, AutoDock-GPU writes the ligand atoms and the co-evaluated
#      flexible-residue atoms (BEGIN_RES/END_RES blocks) inside the SAME
#      DOCKED: MODEL ... ENDMDL block — so a single block-level extraction
#      naturally yields the merged ligand+flexres pose, with no separate
#      stitching step required.
extract_best_pose() {
    local dlg="$1"
    local out_pdbqt="$2"
    local receptor_label="$3"
    local ligand_label="$4"

    python3 - "$dlg" "$out_pdbqt" "$receptor_label" "$ligand_label" <<'PYEOF'
import re
import sys

dlg_path, out_path, receptor_label, ligand_label = sys.argv[1:5]

# Record types that belong in a clean, standalone PDBQT pose file.
KEEP_PREFIXES = ("REMARK", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH",
                  "ATOM", "HETATM", "TORSDOF", "BEGIN_RES", "END_RES", "TER")

with open(dlg_path, "r", errors="replace") as fh:
    text = fh.read()

# ── Split into per-run DOCKED: MODEL ... ENDMDL blocks ──────────────────────
blocks = re.findall(
    r"^DOCKED: MODEL\s+(\d+)\n(.*?)^DOCKED: ENDMDL\n",
    text, flags=re.MULTILINE | re.DOTALL
)
if not blocks:
    print("ERROR: no DOCKED: MODEL ... ENDMDL blocks found in DLG", file=sys.stderr)
    sys.exit(1)

energy_re = re.compile(r"Estimated Free Energy of Binding\s*=\s*([+-]?\d+\.\d+)")

best_run, best_energy, best_block = None, None, None
for run_str, block in blocks:
    m = energy_re.search(block)
    if not m:
        continue
    energy = float(m.group(1))
    if best_energy is None or energy < best_energy:
        best_energy, best_run, best_block = energy, int(run_str), block

if best_block is None:
    print("ERROR: could not parse binding energies from any DOCKED block", file=sys.stderr)
    sys.exit(1)

# ── Strip the "DOCKED: " prefix; keep only recognised PDBQT record types ────
pose_lines = []
for line in best_block.splitlines():
    if not line.startswith("DOCKED: "):
        continue
    content = line[len("DOCKED: "):]
    tokens = content.split()
    tag = tokens[0] if tokens else ""
    if tag in KEEP_PREFIXES:
        pose_lines.append(content.rstrip())

if not pose_lines:
    print("ERROR: best-energy block contained no recognised PDBQT records", file=sys.stderr)
    sys.exit(1)

# ── QC: cross-check ligand vs. flexres atom counts against the DLG header ───
m_lig  = re.search(r"Number of ligand atoms:\s+(\d+)", text)
m_flex = re.search(r"Number of flexres atoms:\s+(\d+)", text)
n_lig_atoms_expected  = int(m_lig.group(1))  if m_lig  else None
n_flex_atoms_expected = int(m_flex.group(1)) if m_flex else None

first_begin_res_idx = next(
    (i for i, l in enumerate(pose_lines) if l.startswith("BEGIN_RES")), None
)

def count_atoms(lines):
    return sum(1 for l in lines if l.startswith(("ATOM", "HETATM")))

if first_begin_res_idx is not None:
    n_lig_atoms_found  = count_atoms(pose_lines[:first_begin_res_idx])
    n_flex_atoms_found = count_atoms(pose_lines[first_begin_res_idx:])
else:
    n_lig_atoms_found  = count_atoms(pose_lines)
    n_flex_atoms_found = 0

qc_notes = []
if n_lig_atoms_expected is not None and n_lig_atoms_found != n_lig_atoms_expected:
    qc_notes.append(f"ligand atom count mismatch (expected {n_lig_atoms_expected}, found {n_lig_atoms_found})")
if n_flex_atoms_expected is not None and n_flex_atoms_found != n_flex_atoms_expected:
    qc_notes.append(f"flexres atom count mismatch (expected {n_flex_atoms_expected}, found {n_flex_atoms_found})")

n_flex_res   = sum(1 for l in pose_lines if l.startswith("BEGIN_RES"))
flex_res_ids = [l.split(None, 1)[1].strip() for l in pose_lines if l.startswith("BEGIN_RES")]

header = [
    "REMARK  Best pose extracted from AutoDock-GPU DLG (gpu_rigid_colab_docking.sh)",
    f"REMARK  Receptor                    : {receptor_label}",
    f"REMARK  Ligand                      : {ligand_label}",
    f"REMARK  Run / Estimated Free Energy : Run {best_run}  |  {best_energy:.2f} kcal/mol",
    f"REMARK  Ligand atoms                : {n_lig_atoms_found}",
]
if n_flex_res:
    header.append(f"REMARK  Flexible residues ({n_flex_res})       : " + ", ".join(flex_res_ids))
    header.append(f"REMARK  Flexres atoms               : {n_flex_atoms_found}")
for note in qc_notes:
    header.append(f"REMARK  QC WARNING                  : {note}")

with open(out_path, "w") as fh:
    fh.write("\n".join(header) + "\n")
    fh.write("\n".join(pose_lines) + "\n")

# Machine-readable summary line consumed by the calling bash function
print(f"RUN={best_run} ENERGY={best_energy:.2f} LIG_ATOMS={n_lig_atoms_found} "
      f"FLEX_RES={n_flex_res} FLEX_ATOMS={n_flex_atoms_found} "
      f"QC_WARNINGS={len(qc_notes)}")
sys.exit(0)
PYEOF
}


# ═══════════════════════════════════════════════════════════
# SECTION 6 — PRE-FLIGHT CHECKS
# ═══════════════════════════════════════════════════════════
echo "── PRE-FLIGHT CHECKS ─────────────────────────────────"

# [R2] GPU / CUDA availability
if ! command -v nvidia-smi &>/dev/null || ! nvidia-smi &>/dev/null; then
    echo "ERROR: No CUDA GPU detected."
    echo "       In Colab: Runtime → Change runtime type → GPU"
    exit 1
fi
echo "✔  GPU     : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"

# [R2] Auto-select the best available AutoDock-GPU binary
AUTODOCK_BIN=""
for wi in 256wi 128wi 64wi; do
    candidate="$BIN_DIR/autodock_gpu_${wi}"
    if [[ -x "$candidate" ]]; then
        AUTODOCK_BIN="$candidate"
        echo "✔  Binary  : autodock_gpu_${wi}"
        break
    fi
done
if [[ -z "$AUTODOCK_BIN" ]]; then
    echo "ERROR: No AutoDock-GPU binary found in $BIN_DIR"
    echo "       Expected: autodock_gpu_256wi / 128wi / 64wi"
    exit 1
fi

# [R3] Log AutoDock-GPU version for publication reproducibility
echo "✔  Version :"
"$AUTODOCK_BIN" --version 2>&1 | head -3 | sed 's/^/   /'

# [R4] python3 is used only for best-pose extraction on flexres receptors;
#      its absence does not block docking, only that one post-processing step.
PYTHON3_AVAILABLE=true
if ! command -v python3 &>/dev/null; then
    PYTHON3_AVAILABLE=false
    echo "⚠  python3   : not found — best-pose extraction will be skipped for"
    echo "              any --flexres receptors (docking itself is unaffected)."
else
    echo "✔  python3  : $(python3 --version 2>&1)"
fi

# [R2] Required directory existence
for dir_label in "RECEPTOR_ROOT_DIR:$RECEPTOR_ROOT_DIR" "LIGAND_DIR:$LIGAND_DIR"; do
    label="${dir_label%%:*}"
    path="${dir_label##*:}"
    if [[ ! -d "$path" ]]; then
        echo "ERROR: Directory not found — $label = $path"
        exit 1
    fi
done

# [R2] Disk space pre-flight
available_gb=$(df -BG "$OUTPUT_ROOT" | awk 'NR==2 { gsub("G",""); print $4 }')
if (( available_gb < MIN_DISK_GB )); then
    echo "⚠  WARNING : Only ${available_gb}GB free (minimum: ${MIN_DISK_GB}GB)"
    echo "   Risk of incomplete run. Free space before continuing."
else
    echo "✔  Disk    : ${available_gb}GB available"
fi

# [R2] Ligand PDBQT discovery
mapfile -t ligand_files < <(find "$LIGAND_DIR" -maxdepth 1 -name "*.pdbqt" 2>/dev/null | sort)
if [[ ${#ligand_files[@]} -eq 0 ]]; then
    echo "ERROR: No .pdbqt ligand files found in $LIGAND_DIR"
    exit 1
fi
echo "✔  Ligands : ${#ligand_files[@]} .pdbqt file(s) found"

# [R2] Receptor subdirectory discovery
mapfile -t receptor_dirs < <(find "$RECEPTOR_ROOT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
if [[ ${#receptor_dirs[@]} -eq 0 ]]; then
    echo "ERROR: No receptor subdirectories found in $RECEPTOR_ROOT_DIR"
    echo "       Expected structure: rigid/<receptor_name>/<receptor>.maps.fld"
    exit 1
fi
echo "✔  Receptors : ${#receptor_dirs[@]} subfolder(s) found"
echo


# ═══════════════════════════════════════════════════════════
# SECTION 7 — LIGAND PDBQT PRE-VALIDATION
#   All ligands are validated once up-front; only valid ones
#   enter the docking loop, saving wasted GPU time.
# ═══════════════════════════════════════════════════════════
echo "── LIGAND PDBQT VALIDATION ───────────────────────────"
valid_ligands=()
invalid_count=0

for ligand in "${ligand_files[@]}"; do
    ligand_name=$(basename "$ligand")
    if validate_pdbqt "$ligand"; then
        valid_ligands+=("$ligand")
    else
        echo "   ⚠  Invalid PDBQT — skipping: $ligand_name"
        echo "N/A,$ligand_name,invalid_pdbqt,N/A,N/A,PDBQT format check failed,N/A,N/A" >> "$REPORT_CSV"
        invalid_count=$(( invalid_count + 1 ))
    fi
done

if [[ ${#valid_ligands[@]} -eq 0 ]]; then
    echo "ERROR: No valid ligand PDBQT files to dock. Exiting."
    exit 1
fi
echo "✔  Valid ligands : ${#valid_ligands[@]} / ${#ligand_files[@]}  (${invalid_count} rejected)"
echo


# ═══════════════════════════════════════════════════════════
# SECTION 8 — DOCKING LOOP
# ═══════════════════════════════════════════════════════════
shopt -s nullglob

# Run counters
SUCCESS=0
FAIL=0
SKIP=0
WARN=0
RESUME=0
FLEXRES_RECEPTORS=0
POSE_EXTRACTED=0

echo "── DOCKING RUN ────────────────────────────────────────"

for receptor_dir in "${receptor_dirs[@]}"; do
    receptor_name=$(basename "$receptor_dir")
    echo
    echo "══ Receptor : $receptor_name"
    echo "   Path     : $receptor_dir"

    # ── [R3] Receptor PDBQT presence warning ────────────────
    #      Excludes *_flex.pdbqt explicitly — both files share the receptor_
    #      name prefix and end in .pdbqt, so an unqualified glob could
    #      nondeterministically report the flexres file here instead of the
    #      rigid receptor file (find's return order isn't alphabetical).
    receptor_pdbqt=$(find "$receptor_dir" -maxdepth 1 -name "*.pdbqt" ! -name "*_flex.pdbqt" 2>/dev/null | head -1)
    if [[ -z "$receptor_pdbqt" ]]; then
        echo "   ⚠  No receptor .pdbqt found — downstream visualization"
        echo "      and re-scoring tools may fail for this receptor."
    else
        echo "   ✔  Receptor PDBQT : $(basename "$receptor_pdbqt")"
    fi

    # ── [R4] Universal --flexres auto-detection ──────────────
    #      Reset every iteration so a flexres receptor never leaks its flag
    #      into the next (rigid-only) receptor folder.
    flexres_pdbqt=""
    flexres_flag="no"
    candidate_flexres="$receptor_dir/${receptor_name}_flex.pdbqt"
    if [[ -f "$candidate_flexres" ]]; then
        if validate_flexres_pdbqt "$candidate_flexres"; then
            flexres_pdbqt="$(realpath "$candidate_flexres")"
            flexres_flag="yes"
            n_flex_res_detected=$(grep -c "^BEGIN_RES" "$candidate_flexres")
            echo "   ✔  Flexres PDBQT  : $(basename "$flexres_pdbqt")  ($n_flex_res_detected flexible residue(s) — --flexres will be added)"
            FLEXRES_RECEPTORS=$(( FLEXRES_RECEPTORS + 1 ))
        else
            echo "   ⚠  ${receptor_name}_flex.pdbqt found but failed structural validation"
            echo "      (expected BEGIN_RES/END_RES/ROOT/ATOM records from prepare_flexreceptor4.py)"
            echo "      Falling back to RIGID-ONLY docking for this receptor."
        fi
    else
        echo "   ℹ  No ${receptor_name}_flex.pdbqt found — rigid receptor docking"
    fi

    # ── [R1] Auto-discover .maps.fld (any filename) ──────────
    mapfile -t fld_files < <(find "$receptor_dir" -maxdepth 1 -name "*.maps.fld" 2>/dev/null)

    if [[ ${#fld_files[@]} -eq 0 ]]; then
        echo "   ⚠  Skipping — no *.maps.fld found in:"
        echo "      $receptor_dir"
        echo "$receptor_name,ALL,skipped_no_fld,N/A,N/A,No .maps.fld file found,$flexres_flag,N/A" >> "$REPORT_CSV"
        SKIP=$(( SKIP + 1 ))
        continue
    fi

    if [[ ${#fld_files[@]} -gt 1 ]]; then
        echo "   ⚠  Skipping — multiple *.maps.fld files detected (ambiguous):"
        for f in "${fld_files[@]}"; do echo "      $f"; done
        echo "      Keep only one .maps.fld per receptor folder."
        echo "$receptor_name,ALL,skipped_ambiguous_fld,N/A,N/A,Multiple .maps.fld files,$flexres_flag,N/A" >> "$REPORT_CSV"
        SKIP=$(( SKIP + 1 ))
        continue
    fi

    # Use realpath to guarantee absolute path survives cd into output dir
    maps_fld="$(realpath "${fld_files[0]}")"
    echo "   ✔  FLD file : $(basename "$maps_fld")"

    # ── [R3] FLD cross-validation (all referenced .map files) ─
    if ! validate_fld "$maps_fld"; then
        echo "   ⚠  Skipping $receptor_name — incomplete map set (details above)"
        echo "$receptor_name,ALL,skipped_incomplete_maps,N/A,N/A,Missing .map files referenced in FLD,$flexres_flag,N/A" >> "$REPORT_CSV"
        SKIP=$(( SKIP + 1 ))
        continue
    fi
    echo "   ✔  Map integrity : all referenced .map files verified"

    receptor_out="$OUTPUT_ROOT/$receptor_name"
    mkdir -p "$receptor_out"

    # ── Per-ligand docking ───────────────────────────────────
    for ligand in "${valid_ligands[@]}"; do
        ligand_name=$(basename "$ligand")
        ligand_stem="${ligand_name%.pdbqt}"
        expected_dlg="$receptor_out/${ligand_stem}.dlg"

        # [R3] Resume: skip pairs already successfully docked
        if [[ -f "$expected_dlg" ]] && grep -q "DOCKED:" "$expected_dlg" 2>/dev/null; then
            echo "   ⏭   Resumed : $ligand_name  (valid DLG already exists)"
            RESUME=$(( RESUME + 1 ))
            SUCCESS=$(( SUCCESS + 1 ))
            continue
        fi

        echo "   ▶   Docking : $ligand_name"
        cp "$ligand" "$receptor_out/"

        # [R4] Universal --flexres injection — only added when a validated
        #      receptor_name_flex.pdbqt was detected for this receptor;
        #      rigid-only receptors get the exact same command as before.
        extra_args=()
        if [[ -n "$flexres_pdbqt" ]]; then
            extra_args+=( --flexres "$flexres_pdbqt" )
        fi

        # [R3] Run docking with per-ligand timeout.
        #      AutoDock-GPU stdout/stderr is captured to a per-ligand tempfile
        #      so that only the concise run-time summary lines reach the log;
        #      the verbose CUDA setup, generation table, and per-run evaluations
        #      are discarded, keeping log size tractable for large campaigns.
        dock_tmplog=$(mktemp)
        (
            cd "$receptor_out"
            timeout "$TIMEOUT_SECONDS" "$AUTODOCK_BIN" \
                -lfile "$ligand_name" \
                -ffile "$maps_fld" \
                "${extra_args[@]}" \
                -nrun "$NRUN" \
                --heuristics 1 \
                --heurmax "$HEURMAX" \
                --nev "$NEV" \
                --lsmet "$LSMET" \
                --lsit "$LSIT" \
                --autostop "$AUTOSTOP" \
                --asfreq "$ASFREQ" \
                --stopstd "$STOPSTD" \
                --clustering "$CLUSTERING" \
                --gbest "$GBEST" \
                --xmloutput "$XML" \
                --dlgoutput "$DLG" \
                -p "$P" \
                -g "$G"
        ) > "$dock_tmplog" 2>&1
        dock_exit=$?

        # Emit only the two concise summary lines produced by AutoDock-GPU
        grep -E "^Run time of entire job set|^All jobs" "$dock_tmplog" | sed 's/^/   /'
        rm -f "$dock_tmplog"

        # Clean up copied ligand from output directory
        rm -f "$receptor_out/$ligand_name"

        # ── [R2+R3] Output verification & error classification ─

        # Timeout
        if [[ $dock_exit -eq 124 ]]; then
            echo "   ⏱   TIMEOUT : $ligand_name exceeded ${TIMEOUT_SECONDS}s"
            echo "$receptor_name,$ligand_name,timeout,N/A,$expected_dlg,Exceeded ${TIMEOUT_SECONDS}s limit,$flexres_flag,N/A" >> "$REPORT_CSV"
            FAIL=$(( FAIL + 1 ))
            continue
        fi

        # Non-zero exit from AutoDock-GPU itself
        if [[ $dock_exit -ne 0 ]]; then
            echo "   ✖   FAILED  : $ligand_name  (exit code $dock_exit)"
            echo "$receptor_name,$ligand_name,failed,N/A,$expected_dlg,Exit code $dock_exit,$flexres_flag,N/A" >> "$REPORT_CSV"
            FAIL=$(( FAIL + 1 ))
            continue
        fi

        # DLG not produced or empty
        if [[ ! -f "$expected_dlg" ]] || [[ ! -s "$expected_dlg" ]]; then
            echo "   ✖   FAILED  : $ligand_name — DLG not produced or is empty"
            echo "$receptor_name,$ligand_name,no_dlg,N/A,$expected_dlg,DLG missing or empty,$flexres_flag,N/A" >> "$REPORT_CSV"
            FAIL=$(( FAIL + 1 ))
            continue
        fi

        # DLG present but contains no pose records
        if ! grep -q "DOCKED:" "$expected_dlg"; then
            echo "   ⚠   WARNING : $ligand_name — DLG written but contains no DOCKED poses"
            echo "$receptor_name,$ligand_name,no_poses,N/A,$expected_dlg,DLG exists but no DOCKED entries,$flexres_flag,N/A" >> "$REPORT_CSV"
            WARN=$(( WARN + 1 ))
            continue
        fi

        # ── [R3] Energy extraction + sanity check ─────────────
        best_energy=$(get_best_energy "$expected_dlg")
        check_energy_sanity "$best_energy" "$ligand_name"

        # ── [R4] Universal flexres-aware best-pose extraction ──
        #      Only attempted when this receptor had a validated flexres file
        #      (rigid-only receptors are completely unaffected). Pulls the
        #      lowest-energy Run out of the DLG and writes the ligand pose
        #      merged with its co-evaluated flexible-residue pose.
        best_pose_pdbqt="N/A"
        if [[ -n "$flexres_pdbqt" ]]; then
            if [[ "$PYTHON3_AVAILABLE" == true ]]; then
                candidate_best_pose="$receptor_out/${ligand_stem}_best_pose.pdbqt"
                extract_summary="$(extract_best_pose "$expected_dlg" "$candidate_best_pose" "$receptor_name" "$ligand_stem" 2>&1)"
                extract_exit=$?
                if [[ $extract_exit -eq 0 ]] && [[ -s "$candidate_best_pose" ]]; then
                    best_pose_pdbqt="$candidate_best_pose"
                    POSE_EXTRACTED=$(( POSE_EXTRACTED + 1 ))
                    echo "   ✔   Best pose (ligand+flexres) → $(basename "$candidate_best_pose")"
                    echo "       $extract_summary"

                    # QC cross-check: extracted-pose energy must match best_energy above —
                    # both use the identical "lowest Estimated Free Energy of Binding" rule.
                    py_energy=$(echo "$extract_summary" | grep -oP 'ENERGY=\K[+-]?[0-9.]+')
                    if [[ -n "$py_energy" ]] && [[ "$best_energy" != "N/A" ]]; then
                        if ! awk "BEGIN { d=($py_energy)-($best_energy); if (d<0) d=-d; exit !(d < 0.01) }"; then
                            echo "   ⚠  QC: extracted-pose energy ($py_energy) differs from reported best_energy ($best_energy) — investigate"
                        fi
                    fi
                else
                    echo "   ⚠   Best-pose extraction failed for $ligand_name (DLG remains the authoritative result):"
                    echo "$extract_summary" | sed 's/^/       /'
                fi
            else
                echo "   ⚠   Skipping best-pose extraction for $ligand_name — python3 unavailable"
            fi
        fi

        echo "   ✔   Done    : $ligand_name  |  Best energy: ${best_energy} kcal/mol"
        echo "$receptor_name,$ligand_name,success,$best_energy,$expected_dlg,,$flexres_flag,$best_pose_pdbqt" >> "$REPORT_CSV"
        SUCCESS=$(( SUCCESS + 1 ))

    done  # end ligand loop
done  # end receptor loop


# ═══════════════════════════════════════════════════════════
# SECTION 9 — FINAL SUMMARY
# ═══════════════════════════════════════════════════════════
echo
echo "══════════════════════════════════════════════════════"
echo "  DOCKING COMPLETE — $(date)"
echo "══════════════════════════════════════════════════════"
echo "  Receptors skipped       : $SKIP"
echo "  Ligands succeeded       : $SUCCESS  (incl. $RESUME resumed)"
echo "  Ligands failed/timeout  : $FAIL"
echo "  Warnings (no poses)     : $WARN"
echo "  Flexres receptors       : $FLEXRES_RECEPTORS  (--flexres auto-applied)"
echo "  Best poses extracted    : $POSE_EXTRACTED  (ligand+flexres merged PDBQT)"
echo "──────────────────────────────────────────────────────"
echo "  Results directory : $OUTPUT_ROOT"
echo "  Full log          : $LOG"
echo "  CSV report        : $REPORT_CSV"
echo "══════════════════════════════════════════════════════"
