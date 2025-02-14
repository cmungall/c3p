"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: CHEBI:36973 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound where monosaccharide units are joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for monosaccharide rings
    ring_atoms = mol.GetRingInfo().AtomRings()
    monosaccharide_rings = []
    for ring in ring_atoms:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if all(atom.GetAtomicNum() in [6, 8] for atom in ring_atoms):  # Carbon and oxygen atoms
            monosaccharide_rings.append(ring)
    
    if not monosaccharide_rings:
        return False, "No monosaccharide rings found"
    
    # Look for glycosidic linkages (acetal bonds between rings)
    for i, ring1 in enumerate(monosaccharide_rings):
        for j, ring2 in enumerate(monosaccharide_rings):
            if i >= j:
                continue
            for atom1 in ring1:
                for atom2 in ring2:
                    if mol.GetBondBetweenAtoms(atom1, atom2):
                        bond = mol.GetBondBetweenAtoms(atom1, atom2)
                        if bond.GetBondType() == Chem.BondType.SINGLE:
                            # Potential glycosidic linkage
                            continue
    
    # Count rotatable bonds (oligosaccharides typically have few rotatable bonds)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 10:
        return False, "Too many rotatable bonds for an oligosaccharide"
    
    return True, "Contains monosaccharide rings connected by glycosidic linkages"