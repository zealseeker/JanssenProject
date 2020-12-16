import warnings

from rdkit import Chem
import pandas as pd
from chembl_structure_pipeline import standardizer

def read_standardise(smi) -> Chem.rdchem.Mol:
    null_mol = Chem.MolFromSmiles('')
    try:
        m = Chem.MolFromSmiles(smi)
        assert m is not None
    except:
        warnings.warn('SMILES: {} is invalid'.format(smi))
        return null_mol
    m = standardizer.standardize_mol(m)
    m = standardizer.get_parent_mol(m)[0]
    return m

def standardise_dataframe(df):
    mols = df['Compound SMILES'].apply(read_standardise)
    df['std_canonical_smiles'] = mols.apply(Chem.MolToSmiles)
    df['std_inchi_key'] = mols.apply(Chem.MolToInchiKey)
    return df

if __name__ == "__main__":
    df = pd.read_csv('../LINCS/input/janssen_compounds.csv')
    df = standardise_dataframe(df)
    df.to_csv('janssen_compounds_std.csv')
