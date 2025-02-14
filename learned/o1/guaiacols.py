"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    
    A guaiacol is defined as any phenol carrying an additional methoxy substituent at the ortho-position.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Flag to indicate if guaiacol substructure is found
    guaiacol_found = False
    
    # Iterate over each ring in the molecule
    for ring in atom_rings:
        ring_set = set(ring)
        # Create a set for faster lookup
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Iterate over atoms in the ring to find phenolic hydroxyl groups
        for atom in ring_atoms:
            # Check if atom is an aromatic carbon
            if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
                # Check if atom is bonded to a hydroxyl group
                is_pheno_hydroxyl = False
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        # Oxygen atom with only one bond (hydroxyl)
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            is_pheno_hydroxyl = True
                            hydroxyl_idx = neighbor.GetIdx()
                            break
                if is_pheno_hydroxyl:
                    # Look for ortho methoxy group
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in ring_set and neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                            # Neighboring aromatic carbon in the ring
                            for sub_neighbor in neighbor.GetNeighbors():
                                if sub_neighbor.GetAtomicNum() == 8:
                                    # Oxygen atom
                                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), sub_neighbor.GetIdx())
                                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                        if sub_neighbor.GetDegree() == 2:
                                            # Oxygen connected to another atom (potential methoxy)
                                            for connected_atom in sub_neighbor.GetNeighbors():
                                                if connected_atom.GetAtomicNum() == 6 and connected_atom.GetIdx() != neighbor.GetIdx():
                                                    if connected_atom.GetDegree() == 1:
                                                        # Methyl group (carbon with one bond)
                                                        guaiacol_found = True
                                                        return True, "Molecule is a guaiacol (phenol with ortho methoxy group)"
    
    if not guaiacol_found:
        return False, "Molecule does not contain guaiacol substructure"