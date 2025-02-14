"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is defined as 'Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.'
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    found = False  # Flag to indicate if tetrahydrofuranone ring is found

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Iterate over all rings in the molecule
    for ring in rings:
        # Consider only five-membered rings
        if len(ring) != 5:
            continue
        
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Check if ring contains exactly one oxygen atom and four carbon atoms
        oxygen_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        carbon_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 6]

        if len(oxygen_atoms) != 1 or len(carbon_atoms) != 4:
            continue  # Not a ring with exactly one oxygen and four carbons
        
        # Check that all bonds in the ring are single (saturated ring)
        is_saturated = True
        for i in range(len(ring)):
            atom1 = mol.GetAtomWithIdx(ring[i])
            atom2 = mol.GetAtomWithIdx(ring[(i+1)%len(ring)])
            bond = mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx())
            if bond.GetBondType() != Chem.BondType.SINGLE:
                is_saturated = False
                break
        if not is_saturated:
            continue  # Ring is not saturated
        
        # Now check for oxo substituent (C=O) on any ring carbon
        for carbon in carbon_atoms:
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring:
                    bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        found = True
                        break
            if found:
                break
        if found:
            break

    if found:
        return True, "Contains an oxolane ring with an oxo substituent (tetrahydrofuranone)"
    else:
        return False, "Does not contain a tetrahydrofuranone structure"