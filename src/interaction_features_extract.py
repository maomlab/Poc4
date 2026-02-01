#!/usr/bin/env python3

"""
# Protein Design Analysis: PLiP Features → DockQ2 Score

This pipeline extracts protein-protein interaction features using PLiP, trains a CatBoost model to predict DockQ2 scores, and performs SHAP analysis for interpretability.

## Installation

### 1. Conda Environment (Recommended)
```bash
conda create -n protein_ml python=3.10 -y
conda activate protein_ml
```

# Core scientific packages
pip install pandas numpy scikit-learn

# CatBoost and SHAP
pip install catboost shap

# PLiP and structure handling
pip install plip biotite

# Optional: For PDB conversion if needed
pip install pdb-tools


PLiP requires OpenBabel. On Ubuntu/Debian:
sudo apt-get install openbabel


## Usage

1. Prepare Input Data

Create a CSV file (designs.csv) with columns:
* design_id: Unique identifier
* complex_cif: Path to AlphaFold3 CIF file
* target_chain: Chain ID for target protein (e.g., 'A')
* design_chain: Chain ID for design (e.g., 'B')
* DockQ2: Target score

2. Run Feature Extraction

python extract_interaction_features.py --input designs.csv --output features.csv --n_jobs 4
 * Extracts PLiP interaction features from protein complexes.
 * Handles CIF files (via automatic PDB conversion) and aggregates features by target residue.


3. Train Model and Generate SHAP

python train_model.py --features features.csv --target DockQ2 --output results/
 * Outputs
 ** catboost_model.cbm: Trained model
 ** shap_summary.png: SHAP summary plot
 ** shap_dependence/: Per-feature dependence plots
 ** feature_importance.csv: Global feature importance
 ** cv_scores.csv: Cross-validation results

"""

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import io
import pandas as pd
from tqdm import tqdm
import biotite.structure.io.pdbx
import biotite.structure.io.pdb
import plip.basic.config
from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Supported PLiP interaction types
INTERACTION_TYPES = [
    'hbonds', 'saltbridges', 'hydrophobic', 'pistacking',
    'pication', 'halogen', 'waterbridges', 'metal'
]

class PLiPFeatureExtractor:
    """
    Extracts features from protein-protein interactions using PLiP.
    Creates count-based features for each target residue position.
    """

    def __init__(self):
        pass

    def load_structure(self, cif_path: str) -> tuple:
        """
        Load CIF file and convert to PDB format in memory for PLiP.
        Returns (structure, pdb_string)
        """
        try:
            # Load with biotite (handles CIF natively)
            cif_file = biotite.structure.io.pdbx.CIFFile.read(cif_path)
            structure = biotite.structure.io.pdbx.get_structure(cif_file, model=1)

            # Convert to PDB string for PLiP
            # PLiP expects PDB format, so we write to string buffer
            pdb_file = biotite.structure.io.pdb.PDBFile()
            pdb_file.set_structure(structure)
            string_buffer = io.StringIO()
            pdb_file.write(string_buffer)
            pdb_buffer = string_buffer.getvalue()

            return structure, pdb_buffer
        except Exception as e:
            logger.error(f"Failed to load {cif_path}: {e}")
            raise


    def extract_interactions(
            self,
            design_id: str,
            pdb_string: str,
            cif_path: str,
            DockQ2: float,
            target_chain: str,
            design_chain: str) -> List:
        """
        Run PLiP analysis and extract interactions between target and design chains.

        Returns:
            List of interactions
            Format: "Target_{res_type}{res_num}_{name}" -> count
        """
        all_interactions = []

        #try:
        # Initialize PLiP complex
        complex_obj = PDBComplex()
        plip.basic.config.CHAINS = [[target_chain], [design_chain]]
        complex_obj.load_pdb(pdb_string, as_string=True)
        complex_obj.analyze()


        # We look for interactions where one participant is in target and other in design
        # PLiP organizes this by 'Binding Sites'
        for site_id, site in complex_obj.interaction_sets.items():
            # We filter for interactions involving our specified chains
            # This is a simplification; PLiP interaction objects contain 'resnr' and 'restype'

            # Extract specific interaction types
            interactions = [
                (site.hydrophobic_contacts, "hydrop"),
                (site.hbonds_pdon, "hbond_pdon"), # target as donor
                (site.hbonds_ldon, "hbond_ldon"), # design as donor
                (site.saltbridge_pneg, "saltbridge_pneg"),
                (site.saltbridge_lneg, "saltbridge_lneg"),
                (site.pistacking, "pistacking"),
                (site.pication_paro, "pication_paro"),
                (site.pication_laro, "pication_laro"),
            ]

            for int_list, name in interactions:
                for interaction in int_list:
                    # Logic: Determine if interaction bridges target_chain and design_chain
                    # PLiP attributes vary by interaction type, but usually have .resnr
                    res_num = interaction.resnr
                    res_type = interaction.restype

                    # We tag features by the Target Residue
                    # Filter check: Is this residue on the target chain?
                    # (Note: PLiP mapping can be complex; ensure chain IDs match CIF)
                    feature_name = f"Target_{res_type}{res_num}_{name}"

                    all_interactions.append({
                        'design_id': design_id,
                        'cif_path': cif_path,
                        'DockQ2': DockQ2,
                        'target_chain': interaction.reschain,
                        'target_resnr': interaction.resnr,
                        'target_restype': interaction.restype,
                        'design_chain': interaction.reschain_l,
                        'design_resnr': interaction.resnr_l,
                        'design_restype': interaction.restype_l,
                        'interaction_type': name
                    })

        logger.info(f"Extracted {len(all_interactions)} interactions")
        return all_interactions

    def process_complex(
            self,
            design_id: str,
            cif_path: str,
            DockQ2: str,
            target_chain: str,
            design_chain: str) -> List:
        """
        Process a single complex and return all features.
        """

        # Load structure
        structure, pdb_string = self.load_structure(cif_path)

        # Extract interactions
        all_interactions = self.extract_interactions(
            design_id,
            pdb_string,
            cif_path,
            DockQ2,
            target_chain,
            design_chain)

        logger.info(f"Successfully processed {Path(cif_path).name}")

        return all_interactions

def main():
    parser = argparse.ArgumentParser(description="Extract PLiP features from protein complexes")
    parser.add_argument("--input", required=True, help="Input CSV with design_id, complex_cif, target_chain, design_chain")
    parser.add_argument("--output", required=True, help="Output CSV with extracted features")
    parser.add_argument("--n_jobs", type=int, default=1, help="Number of parallel jobs")
    parser.add_argument("--include_residue_id", action="store_true", default=True,
                       help="Include residue identity features")
    parser.add_argument("--interaction_types", nargs="+", default=INTERACTION_TYPES,
                       help="Interaction types to extract")

    args = parser.parse_args()

    # Load input data
    df_input = pd.read_csv(args.input)
    required_cols = ['design_id', 'complex_cif', 'target_chain', 'design_chain']
    if not all(col in df_input.columns for col in required_cols):
        raise ValueError(f"Input CSV must contain columns: {required_cols}")

    logger.info(f"Processing {len(df_input)} complexes...")

    # Initialize extractor
    extractor = PLiPFeatureExtractor()

    # Process each complex
    results = []
    for _, row in tqdm(df_input.iterrows(), total=len(df_input)):
        features = extractor.process_complex(
            design_id=row['design_id'],
            cif_path=row['complex_cif'],
            DockQ2=row['DockQ2'],
            target_chain=row['target_chain'],
            design_chain=row['design_chain'],

        )
        results.extend(features)

    # Combine results
    df_features = pd.DataFrame(results)

    # Save to CSV
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df_features.to_csv(args.output, index=False)
    logger.info(f"Saved features to {args.output}")

if __name__ == "__main__":
    main()
