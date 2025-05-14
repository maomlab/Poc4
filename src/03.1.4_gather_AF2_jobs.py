
import pathlib
import pickle
import glob
import numpy as np
import pandas as pd
import tqdm

def gather_AF2_jobs(
    job_path):
    rows = []
    pickle_files = glob.glob(f"{job_path}/output/*/input/result_model_*_pred_0.pkl")
    for af2_results_fname in tqdm.tqdm(pickle_files):
        with open(af2_results_fname, "br") as af2_results_file:
            af2_results = pickle.load(af2_results_file)
        rows.append({
            "af2_results_fname" : af2_results_fname,
            "structure_id" : pathlib.Path(af2_results_fname).parent.parent.name,
            "af2_plddt": af2_results["ranking_confidence"],
            "rank" : int(pathlib.Path(af2_results_fname).name.replace("result_model_", "").replace("_pred_0.pkl", ""))
        })
    af2_scores = pd.DataFrame(rows)
    return af2_scores

af2_scores = gather_AF2_jobs(
    job_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0")

af2_scores.to_csv(
    "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/af2_scores.tsv",
    index = False,
    sep = "\t")

    
    
