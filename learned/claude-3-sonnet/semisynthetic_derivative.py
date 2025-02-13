"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: CHEBI:48927 semisynthetic derivative

A semisynthetic derivative is any organic molecular entity derived from a natural product by partial chemical synthesis.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for known natural product scaffolds or building blocks
    murcko_scaffold = Chem.MurckoScaffold.MurckoScaffoldSmilesFromMolBlock(Chem.MolToMolBlock(mol))
    known_scaffolds = ['c1ccc2ccccc2c1', 'c1ccccc1', 'C1CCCCC1', 'C1CCCCCC1', 'C1CCCCCCC1', 'C1CCCCCCCC1', 'C12CCCCC1CCCCC2',
                       'c1cc2ccccc2cc1', 'c1ccc2ccccn2c1', 'c1ccc2ncccn2c1', 'C12C=CC=CC1=CC=CC=C2', 'C12C3C=CC=CC3=CC=CC=C12',
                       'C1CCC2CCCCC2C1', 'C12CC3CC(CC(C3)C1)C2', 'C12C3C4C5C2C1C6C7C3C4C5C67', 'C12C3C4C5C6C7C8C3C41C28C97C56C7']

    if murcko_scaffold in known_scaffolds:
        # Check for additional modifications (e.g. added functional groups, substituents)
        mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
        scaffold_mol = Chem.MolFromSmiles(murcko_scaffold)
        scaffold_weight = rdMolDescriptors.CalcExactMolWt(scaffold_mol)
        if mol_weight > scaffold_weight + 50:  # Arbitrary weight difference threshold
            return True, f"Contains {murcko_scaffold} scaffold and additional modifications"
        else:
            return False, f"Contains {murcko_scaffold} scaffold but no significant modifications"
    else:
        return False, "Does not contain a known natural product scaffold"