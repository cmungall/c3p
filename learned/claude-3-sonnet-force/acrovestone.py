"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: CHEBI:72099 Acrovestone
A polyphenol that is isolated from Acronychia pedunculata and exhibits moderate antioxidant and antityrosinase activities.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone based on its SMILES string.
    Acrovestones are polyphenols with moderate antioxidant and antityrosinase activities, isolated from Acronychia pedunculata.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Acrovestones are polyphenolic compounds, so should contain multiple phenol groups
    phenol_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
    n_phenol_groups = len(mol.GetSubstructMatches(phenol_pattern))
    if n_phenol_groups < 2:
        return False, "Too few phenol groups for acrovestone"

    # Acrovestones often have glycosidic substituents, so look for sugar groups
    has_sugar_groups = any(atom.GetSmarts() in ['OC', 'OCC'] for atom in mol.GetAtoms())
    
    # Calculate common molecular descriptors and properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    tpsa = Descriptors.TPSA(mol)
    logp = Descriptors.MolLogP(mol)

    # Use a combination of structural features and molecular properties to classify
    if (
        n_phenol_groups >= 2
        and 300 < mol_wt < 800
        and n_rotatable_bonds > 5
        and n_rings >= 3
        and n_aromatic_rings >= 2
        and 100 < tpsa < 200
        and 2 < logp < 6
    ):
        if has_sugar_groups:
            return True, "Polyphenolic structure with glycosidic substituents, consistent with acrovestone"
        else:
            return True, "Polyphenolic structure, consistent with acrovestone"
    else:
        return False, "Structural features and molecular properties not typical of acrovestones"