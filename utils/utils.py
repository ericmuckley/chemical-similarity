import os
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem import AllChem, Draw
from sklearn.decomposition import PCA


def get_pca(smiles_list: list):
    """From a list of SMILES strings, get thir 2D PCA representation"""
    # get array of all the molecular fingerprints
    for si, s in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(s)
        fp = Chem.RDKFingerprint(mol)
        bitlist = list(fp.ToBitString())
        if si == 0:
            fp_array = np.empty((len(smiles_list), len(bitlist)))
        fp_array[si] = bitlist

    # perform PCA
    pca_model = PCA(n_components=2)
    pca = pca_model.fit_transform(fp_array)
    # normlize PCA componenets
    pca = (pca - pca.min(0)) / pca.ptp(0)
    return pca.tolist()


def smiles_to_svg(smiles: str):
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
