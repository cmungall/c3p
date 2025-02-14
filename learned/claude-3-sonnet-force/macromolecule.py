"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:24116 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a molecule of high relative molecular mass, the structure
    of which essentially comprises the multiple repetition of units derived,
    actually or conceptually, from molecules of low relative molecular mass.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Macromolecules typically have a molecular weight > 1000 Da
    if mol_wt < 1000:
        return False, f"Molecular weight {mol_wt:.2f} Da is too low for a macromolecule"
    
    # Look for polymeric repeat units (e.g., peptides, polysaccharides, polynucleotides)
    repeat_units = ['C(=O)N', 'OCC', 'OC=O', 'NC=O', 'OP(O)(=O)OC']
    for pattern in repeat_units:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            n_matches = len(mol.GetSubstructMatches(patt))
            if n_matches > 10:  # Arbitrary threshold, adjust as needed
                return True, f"Contains {n_matches} repetitive units derived from smaller molecules"
    
    # Check for long carbon chains
    chain_length = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.SINGLE:
            atoms = [bond.GetBeginAtom(), bond.GetEndAtom()]
            if all(atom.GetAtomicNum() == 6 for atom in atoms):
                chain_length += 1
            else:
                chain_length = 0
        else:
            chain_length = 0
        if chain_length > 20:
            return True, "Contains a long carbon chain characteristic of a macromolecule"
    
    return False, "No characteristic macromolecular features found"