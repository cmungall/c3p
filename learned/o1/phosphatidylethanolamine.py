"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    # Glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    
    # Ester linkage at positions 1 and 2
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    
    # Phosphate group attached to carbon 3
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    
    # Ethanolamine attached to phosphate
    ethanolamine_pattern = Chem.MolFromSmarts("OP(=O)(OCCN)O")
    
    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check for two ester linkages
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester group(s), need at least 2"
    
    # Check for phosphate group at position 3
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found attached to glycerol backbone"
    
    # Check for ethanolamine group attached to phosphate
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "No ethanolamine group found attached to phosphate"
    
    # Confirm that fatty acid chains are long (e.g., more than 6 carbons)
    # Find all ester-linked alkyl chains
    fatty_acid_chains = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.ESTER:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            chain = Chem.PathToSubmol(mol, [atom2.GetIdx()])
            num_carbons = sum(1 for atom in chain.GetAtoms() if atom.GetAtomicNum() == 6)
            fatty_acid_chains.append(num_carbons)
    if len(fatty_acid_chains) < 2 or not all(carbons >= 6 for carbons in fatty_acid_chains):
        return False, "Fatty acid chains are missing or too short"
    
    return True, "Molecule is a phosphatidylethanolamine"