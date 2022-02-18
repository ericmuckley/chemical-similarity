import os
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem import AllChem, Draw
from sklearn.decomposition import PCA


def smiles_to_svg(smiles):
    """Convert a smiles string to an SVG string"""
    mol = Chem.MolFromSmiles(smiles)
    Compute2DCoords(mol)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    drawing = drawer.GetDrawingText()
    return drawing


def get_similarities(smiles0: str, smiles_list: list):
    """
    Calculate the Tanimoto similarity between a single smiles
    string (smiles0) and a list of other smiles strings (smiles_list).
    """
    mol0 = Chem.MolFromSmiles(smiles0)
    fp0 = Chem.RDKFingerprint(mol0)
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    fps = [Chem.RDKFingerprint(m) for m in mols]
    similarities = [DataStructs.FingerprintSimilarity(fp0, fp) for fp in fps]
    return similarities






def draw_molecules(smiles_list, source, img_dir):
    """
    Draw molecules from a list of SMILES strings
    and save the iamges to img_dir.
    """
    # remove everything in the image directory
    #for f in os.listdir(img_dir):
    #    os.remove(os.path.join(img_dir, f))

    # create the images
    for si, s in enumerate(smiles_list):
        img_fp = os.path.join(img_dir, f"{source}-{si}.png")
        mol = Chem.MolFromSmiles(s)
        Draw.MolToFile(mol, img_fp)


def compute_ecfp_descriptors(smiles_list: list):
    """Computes ECPF descriptors"""
    keep_idx = []
    descriptors = []
    for i, smiles in enumerate(smiles_list):
        ecfp = _compute_single_ecfp_descriptor(smiles)
        if ecfp is not None:
            keep_idx.append(i)
            descriptors.append(ecfp)
    return np.vstack(descriptors), keep_idx


def _compute_single_ecfp_descriptor(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as E:
        return None
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        return np.array(fp)
    return None


def standardize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    stan_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    return stan_smiles


def run_pca(X, n_components=2):
    """Perform PCA on input array X and return
    results as a pandas dataframe."""
    pca = PCA(n_components=n_components)
    Xpca = pca.fit_transform(X)
    df = pd.DataFrame({f"pca-{i}": Xpca[:, i] for i in range(Xpca.shape[1])})
    return df, pca