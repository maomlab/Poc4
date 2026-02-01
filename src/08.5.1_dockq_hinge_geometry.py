import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import tqdm
import biotite.structure.io.pdbx as pdbx
import biotite.database.rcsb as rcsb
import biotite.structure as structure


def load_structure(cif_fname, verbose = False):
    if verbose:
        print(f"Loading structure from '{cif_fname}'")

    try:
        struct = pdbx.CIFFile.read(cif_fname)
        struct = pdbx.get_structure(struct, model=1, include_bonds=True)
        return struct
    except Exception as e:
        print(f"Failed to load structure {cif_fname}")
        print(f"ERROR: {e}")


def compute_distance(
        struct,
        chain_id_1, res_id_1, atom_name_1,
        chain_id_2, res_id_2, atom_name_2,
        verbose = False):

    try:
        atom1 = struct[
            (struct.chain_id == chain_id_1) &
            (struct.res_id == res_id_1) &
            (struct.atom_name == atom_name_1)]
        atom2 = struct[
            (struct.chain_id == chain_id_2) &
            (struct.res_id == res_id_2) &
            (struct.atom_name == atom_name_2)]
        distance = structure.distance(atom1, atom2)[0]
        return {'distance' : distance}

    except Exception as e:
        print(f"Failed to compute features for structure")
        print(f"ERROR: {e}")
    return None


def main():

    verbose = True

    intermediate_path = Path("intermediate_data/dockq_analysis/hinge_model")
    intermediate_path.mkdir(parents=True, exist_ok=True)
    date_code = "20260121"

    ##############
    # Train Data #
    ##############

    # find all the cif files to compute features for
    cif_fnames = []
    for cif_fnames_glob in [
            "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_*/seq*/*cif",
            "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_*/fold_seq*/*cif"]:
        cif_fnames.extend(
            glob.glob(cif_fnames_glob))

    if verbose:
        print(f"Computing structural features for {len(cif_fnames)} structures")

    features = []
    for cif_fname in tqdm.tqdm(cif_fnames):
        struct = load_structure(
            cif_fname = cif_fname,
            verbose = verbose)

        if "af3_afs_results_20251003" in cif_fname:
            chain_id_1 = "E"
        else:
            chain_id_1 = "D"

        struct_features = compute_distance(
            struct = struct,
            chain_id_1 = chain_id_1, res_id_1 = 83, atom_name_1 = 'CA', # Pup2
            chain_id_2 = 'K',        res_id_2 = 63, atom_name_2 = 'CA', # Poc4
            verbose = verbose)

        if struct_features is not None:
            struct_features["cif_fname"] = cif_fname
            features.append(struct_features)

    features = pd.DataFrame(features)
    features.to_csv(
        intermediate_path / f"design_hinge_features_{date_code}.tsv",
        sep = "\t",
        index = False)



    #################
    ### Test data ###
    #################

        # find all the cif files to compute features for

    cif_fnames = []
    for cif_fnames_glob in [
            "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/afs_outputs_january/folds_*/*/*.cif",
            "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/afs_outputs_january/fold_*/*.cif"]:
        cif_fnames.extend(
            glob.glob(cif_fnames_glob))

    if verbose:
        print(f"Computing structural features for {len(cif_fnames)} structures")

    features = []
    for cif_fname in tqdm.tqdm(cif_fnames):
        struct = load_structure(
            cif_fname = cif_fname,
            verbose = verbose)

        if "af3_afs_results_20251003" in cif_fname:
            chain_id_1 = "E"
        else:
            chain_id_1 = "D"

        struct_features = compute_distance(
            struct = struct,
            chain_id_1 = chain_id_1, res_id_1 = 83, atom_name_1 = 'CA', # Pup2
            chain_id_2 = 'K',        res_id_2 = 63, atom_name_2 = 'CA', # Poc4
            verbose = verbose)

        if struct_features is not None:
            struct_features["cif_fname"] = cif_fname
            features.append(struct_features)

    features = pd.DataFrame(features)
    features.to_csv(
        intermediate_path / f"design_hinge_features_test_{date_code}.tsv",
        sep = "\t",
        index = False)


if __name__ == "__main__":
    main()
