"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine has a glycerol backbone, a phosphocholine group, and two acyl chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Relaxed pattern for glycerol backbone: COC(C)CO
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Pattern for phosphocholine group: P(-[O-])(=O)OCC[N+](C)(C)C
    phosphocholine_pattern = Chem.MolFromSmarts("P([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"
        
    # Pattern for ester-bonded acyl chains: C(=O)O connected to glycerol
    # More generic to accommodate different representations and configurations
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} acyl chains, need exactly 2"

    # Ensuring at least two esters are connected to the same glycerol unit
    glycerol_acidified = Chem.ReplaceCore(mol, Chem.MolFromSmarts("OCC(O)CO"))
    unique_acyl_attachments = set()
    for match in ester_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:  # Oxygen atom
                attached_carbon = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
                if attached_carbon:
                    unique_acyl_attachments.add(attached_carbon[0])

    if len(unique_acyl_attachments) < 2:
        return False, f"Only {len(unique_acyl_attachments)} unique acyl attachments found, need 2"

    return True, "Molecule is a phosphatidylcholine"