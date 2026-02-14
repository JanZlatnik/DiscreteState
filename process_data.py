#!/usr/bin/env python3
import os
import shutil
import argparse

# ================= CONFIGURATION =================
SOURCE_DIR = "DATA"
OUT_CUT = "DATAcut"
OUT_SPARSE = "DATAsparse"

# Files to specifically process (relative path inside DATA)
TARGET_FILES = {
    os.path.normpath("DSState/defect.txt"),
    os.path.normpath("DSState/delta.txt"),
    os.path.normpath("DSState/deltaA.txt"),
    os.path.normpath("DSState/DS_defect.txt"),
    os.path.normpath("DSState/DS_phaseshift.txt"),
    os.path.normpath("DSState/gamma.txt"),
    os.path.normpath("DSState/gammaA.txt"),
    os.path.normpath("DSState/phaseshift.txt"),
    os.path.normpath("DSState/Vde.txt"),
    os.path.normpath("Hilbert/DeltaContinuous.txt"),
    os.path.normpath("Hilbert/DeltaFull.txt"),
    os.path.normpath("Hilbert/DeltaRydberg.txt"),
    os.path.normpath("Hilbert/Vde.txt")
}
# =================================================

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def process_file(src, dest_cut, dest_sparse, emin, emax, step):
    """Streams file: filters by energy for CUT, skips lines for SPARSE."""
    try:
        with open(src, 'r') as f_in, \
             open(dest_cut, 'w') as f_cut, \
             open(dest_sparse, 'w') as f_sparse:
            
            idx = 0
            for line in f_in:
                # Always keep headers
                if line.strip().startswith('#'):
                    f_cut.write(line)
                    f_sparse.write(line)
                    continue

                try:
                    parts = line.split()
                    if not parts: continue
                    energy = float(parts[0])

                    # CUT: write if within range
                    if emin <= energy <= emax:
                        f_cut.write(line)

                    # SPARSE: write if index matches step
                    if idx % step == 0:
                        f_sparse.write(line)
                    
                    idx += 1

                except ValueError:
                    continue 

    except Exception as e:
        print(f"[ERROR] {src}: {e}")

def main():
    # Argument Parsing
    parser = argparse.ArgumentParser(description="Process DATA files for plotting.")
    parser.add_argument("--emin", type=float, default=-15.0, help="Min Energy (eV) [Default: -15.0]")
    parser.add_argument("--emax", type=float, default=15.0, help="Max Energy (eV) [Default: 15.0]")
    parser.add_argument("--step", type=int, default=1000, help="Sparse step (N-th point) [Default: 1000]")
    
    args = parser.parse_args()

    base_cwd = os.getcwd()
    src_full = os.path.join(base_cwd, SOURCE_DIR)

    if not os.path.exists(src_full):
        print(f"Error: '{SOURCE_DIR}' not found in {base_cwd}")
        return

    print(f"--- Configuration ---")
    print(f"Source: {src_full}")
    print(f"Cut Range: {args.emin} to {args.emax} eV")
    print(f"Sparse Step: {args.step}")
    print(f"---------------------")

    # Processing Loop
    for root, dirs, files in os.walk(src_full):
        rel_path = os.path.relpath(root, src_full)
        
        path_cut = os.path.join(base_cwd, OUT_CUT, rel_path)
        path_sparse = os.path.join(base_cwd, OUT_SPARSE, rel_path)
        
        ensure_dir(path_cut)
        ensure_dir(path_sparse)

        for file in files:
            if not file.endswith(".txt"):
                continue

            src_file = os.path.join(root, file)
            dst_cut_file = os.path.join(path_cut, file)
            dst_sparse_file = os.path.join(path_sparse, file)

            file_rel = os.path.normpath(os.path.join(rel_path, file))
            
            # Check if current file is in the target list
            is_target = any(file_rel == t for t in TARGET_FILES)

            if is_target:
                process_file(src_file, dst_cut_file, dst_sparse_file, args.emin, args.emax, args.step)
                print(f"[PROCESSED] {file_rel}")
            else:
                shutil.copy2(src_file, dst_cut_file)
                shutil.copy2(src_file, dst_sparse_file)

    print("Done.")

if __name__ == "__main__":
    main()