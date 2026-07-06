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
echo "receptor,ligand,status,best_energy_kcal_mol,dlg_path,notes" > "$REPORT_CSV"


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
        echo "N/A,$ligand_name,invalid_pdbqt,N/A,N/A,PDBQT format check failed" >> "$REPORT_CSV"
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

echo "── DOCKING RUN ────────────────────────────────────────"

for receptor_dir in "${receptor_dirs[@]}"; do
    receptor_name=$(basename "$receptor_dir")
    echo
    echo "══ Receptor : $receptor_name"
    echo "   Path     : $receptor_dir"

    # ── [R3] Receptor PDBQT presence warning ────────────────
    receptor_pdbqt=$(find "$receptor_dir" -maxdepth 1 -name "*.pdbqt" 2>/dev/null | head -1)
    if [[ -z "$receptor_pdbqt" ]]; then
        echo "   ⚠  No receptor .pdbqt found — downstream visualization"
        echo "      and re-scoring tools may fail for this receptor."
    else
        echo "   ✔  Receptor PDBQT : $(basename "$receptor_pdbqt")"
    fi

    # ── [R1] Auto-discover .maps.fld (any filename) ──────────
    mapfile -t fld_files < <(find "$receptor_dir" -maxdepth 1 -name "*.maps.fld" 2>/dev/null)

    if [[ ${#fld_files[@]} -eq 0 ]]; then
        echo "   ⚠  Skipping — no *.maps.fld found in:"
        echo "      $receptor_dir"
        echo "$receptor_name,ALL,skipped_no_fld,N/A,N/A,No .maps.fld file found" >> "$REPORT_CSV"
        SKIP=$(( SKIP + 1 ))
        continue
    fi

    if [[ ${#fld_files[@]} -gt 1 ]]; then
        echo "   ⚠  Skipping — multiple *.maps.fld files detected (ambiguous):"
        for f in "${fld_files[@]}"; do echo "      $f"; done
        echo "      Keep only one .maps.fld per receptor folder."
        echo "$receptor_name,ALL,skipped_ambiguous_fld,N/A,N/A,Multiple .maps.fld files" >> "$REPORT_CSV"
        SKIP=$(( SKIP + 1 ))
        continue
    fi

    # Use realpath to guarantee absolute path survives cd into output dir
    maps_fld="$(realpath "${fld_files[0]}")"
    echo "   ✔  FLD file : $(basename "$maps_fld")"

    # ── [R3] FLD cross-validation (all referenced .map files) ─
    if ! validate_fld "$maps_fld"; then
        echo "   ⚠  Skipping $receptor_name — incomplete map set (details above)"
        echo "$receptor_name,ALL,skipped_incomplete_maps,N/A,N/A,Missing .map files referenced in FLD" >> "$REPORT_CSV"
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
            echo "$receptor_name,$ligand_name,timeout,N/A,$expected_dlg,Exceeded ${TIMEOUT_SECONDS}s limit" >> "$REPORT_CSV"
            FAIL=$(( FAIL + 1 ))
            continue
        fi

        # Non-zero exit from AutoDock-GPU itself
        if [[ $dock_exit -ne 0 ]]; then
            echo "   ✖   FAILED  : $ligand_name  (exit code $dock_exit)"
            echo "$receptor_name,$ligand_name,failed,N/A,$expected_dlg,Exit code $dock_exit" >> "$REPORT_CSV"
            FAIL=$(( FAIL + 1 ))
            continue
        fi

        # DLG not produced or empty
        if [[ ! -f "$expected_dlg" ]] || [[ ! -s "$expected_dlg" ]]; then
            echo "   ✖   FAILED  : $ligand_name — DLG not produced or is empty"
            echo "$receptor_name,$ligand_name,no_dlg,N/A,$expected_dlg,DLG missing or empty" >> "$REPORT_CSV"
            FAIL=$(( FAIL + 1 ))
            continue
        fi

        # DLG present but contains no pose records
        if ! grep -q "DOCKED:" "$expected_dlg"; then
            echo "   ⚠   WARNING : $ligand_name — DLG written but contains no DOCKED poses"
            echo "$receptor_name,$ligand_name,no_poses,N/A,$expected_dlg,DLG exists but no DOCKED entries" >> "$REPORT_CSV"
            WARN=$(( WARN + 1 ))
            continue
        fi

        # ── [R3] Energy extraction + sanity check ─────────────
        best_energy=$(get_best_energy "$expected_dlg")
        check_energy_sanity "$best_energy" "$ligand_name"

        echo "   ✔   Done    : $ligand_name  |  Best energy: ${best_energy} kcal/mol"
        echo "$receptor_name,$ligand_name,success,$best_energy,$expected_dlg," >> "$REPORT_CSV"
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
echo "──────────────────────────────────────────────────────"
echo "  Results directory : $OUTPUT_ROOT"
echo "  Full log          : $LOG"
echo "  CSV report        : $REPORT_CSV"
echo "══════════════════════════════════════════════════════"
