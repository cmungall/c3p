"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are naturally occurring compounds derived from prostanoic acid (C20),
    containing a cyclopentane ring connected to two aliphatic chains and functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Carbon count is {c_count}, expected around 20 carbons for prostaglandins"

    # Identify 5-membered rings (cyclopentane rings)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    five_membered_rings = [ring for ring in atom_rings if len(ring) == 5]

    if not five_membered_rings:
        return False, "No 5-membered cyclopentane ring found"

    # Assume prostaglandin if any 5-membered ring matches criteria
    for ring in five_membered_rings:
        # Check for exocyclic attachments (side chains)
        side_chains = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring:
                    side_chains += 1
        # Typically, there should be at least two side chains attached to the ring
        if side_chains >= 2:
            # Check for carboxylic acid group (-COOH)
            carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
            if mol.HasSubstructMatch(carboxylic_acid):
                return True, "Contains cyclopentane ring with side chains and carboxylic acid group typical of prostaglandins"
            else:
                return False, "Carboxylic acid group not found"
    return False, "Cyclopentane ring does not have required side chains"

# Examples of usage
if __name__ == "__main__":
    smiles_list = [
        "CCCCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CCCC(O)=O",  # Prostaglandin F2alpha
        "CC(C)OC(=O)CCC\\C=C/C[C@H]1[C@@H](O)C[C@@H](O)[C@@H]1\\C=C\\C(F)(F)COc1ccccc1"  # Tafluprost
    ]
    for smiles in smiles_list:
        result, reason = is_prostaglandin(smiles)
        print(f"SMILES: {smiles}\nIs prostaglandin: {result}\nReason: {reason}\n")