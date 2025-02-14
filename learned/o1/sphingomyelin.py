"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is defined as a phospholipid where the amino group of a sphingoid base is in amide linkage with a fatty acid,
    and the terminal hydroxyl group is esterified to phosphorylcholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for sphingoid base backbone (long-chain amino alcohol)
    # Define a sphingoid base pattern
    sphingoid_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#7]-[#6]-[#6]-[#8]")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base backbone found"
    
    # Check for amide linkage with fatty acid chain
    # Amide group pattern: carbonyl group connected to nitrogen
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide linkage found between amino group and fatty acid"

    # Check for long-chain fatty acid linked via amide bond
    fatty_acid_found = False
    for match in amide_matches:
        carbonyl_c_idx = match[0]
        amide_n_idx = match[2]
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        amide_n = mol.GetAtomWithIdx(amide_n_idx)
        # Check the length of the carbon chain attached to the carbonyl carbon
        chain_length = 0
        visited_atoms = set()
        atoms_to_visit = [nbr.GetIdx() for nbr in carbonyl_c.GetNeighbors() if nbr.GetIdx() != amide_n_idx]
        while atoms_to_visit:
            atom_idx = atoms_to_visit.pop()
            if atom_idx in visited_atoms:
                continue
            visited_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
                # Add neighbors to visit
                atoms_to_visit.extend([nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in visited_atoms and nbr.GetAtomicNum() == 6])
        if chain_length >= 10:
            fatty_acid_found = True
            break
    if not fatty_acid_found:
        return False, "Amide linkage does not involve a long-chain fatty acid"
    
    # Check for phosphorylcholine group
    phosphorylcholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphorylcholine_pattern):
        return False, "No phosphorylcholine group found"
    
    # Verify that the phosphorylcholine group is attached via ester linkage to the sphingoid base
    ester_linkage_pattern = Chem.MolFromSmarts("[O][P](=O)([O-])OCC[N+](C)(C)C")
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    if not ester_matches:
        return False, "Phosphorylcholine group not attached via ester linkage"
    
    # All checks passed
    return True, "Molecule is a sphingomyelin with sphingoid base, amide-linked fatty acid, and phosphorylcholine group"