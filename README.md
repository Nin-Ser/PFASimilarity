# PFASimilarity
This open-source R script implements a similarity-based framework for matching detected PFAS with structurally similar reference standards.
The approach combines molecular fingerprints (ECFP6 and MACCS keys) and fluorine content derived from SMILES representations to compute weighted Tanimoto similarity scores. For each detected PFAS, the most similar reference standard is identified and reported together with its similarity score.
The workflow is automated and processes two user-provided input files (PFAS and standards lists in CSV format containing compound names and SMILES). It performs molecular parsing, fingerprint calculation, similarity scoring, and outputs a matching table.
Invalid or missing SMILES entries are automatically flagged and excluded from the similarity calculation.
This framework was developed to support the selection of reference standards in PFAS semi-quantification workflows. By providing a reproducible measure of structural similarity between detected compounds and available standards, it helps standardize analogue selection and reduce expert-driven subjectivity.
# Input requirements
-	pfas.csv: detected compounds (Name, SMILES) 
-	standards.csv: reference standards (Name, SMILES) 
# Output
-	similarity.csv: best-matching reference standard for each PFAS compound with similarity score 
# Citation
If you use this script, please cite:
Ninon Serre, Randolph Singh, Lise Boulard, Aurore Zalouk-Vergnoux, Yann Aminot, Towards Reliable Semi-Quantification of PFAS in Aquatic Biota: A Similarity-Based Approach and Practical Recommendations, Talanta, 2026, 130056, https://doi.org/10.1016/j.talanta.2026.130056.

