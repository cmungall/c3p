"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic molecular properties
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20:  # Lowered threshold
        return False, "Too small to be a proanthocyanidin"

    # Count key atoms
    num_o = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_o < 4:  # Lowered threshold
        return False, "Insufficient oxygen atoms"

    # Basic flavan unit patterns (more flexible)
    flavan_patterns = [
        # Basic flavan-3-ol core
        "[#6]1-[#6]-[#6]-c2c(O)cc(O)cc2O1",
        # Alternative pattern with different hydroxylation
        "[#6]1-[#6]-[#6]-c2c(O)cc(O)c(O)c2O1",
        # Pattern for gallated units
        "[#6]1-[#6]-[#6]-c2c(O)cc(O)cc2O[#6]1OC(=O)c1cc(O)c(O)c(O)c1"
    ]
    
    total_flavan_matches = 0
    for pattern in flavan_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            matches = len(mol.GetSubstructMatches(patt))
            total_flavan_matches += matches

    if total_flavan_matches < 1:
        return False, "No flavan unit found"

    # Check for characteristic linkage patterns
    linkage_patterns = [
        # 4→8 linkage
        "[#6]1-[#6]-[#6]-c2c(-[#6]1)cc(O)c1c2O[#6][#6][#6]c2cc(O)cc(O)c21",
        # 4→6 linkage
        "[#6]1-[#6]-[#6]-c2c(-[#6]1)cc(O)c(-[#6]1-[#6]-[#6]-c3c(-[#6]1)cc(O)cc3O)c2O",
        # 2→O→7 linkage (A-type)
        "[#6]1O[#6]2Oc3cc(O)cc(O)c3[#6][#6]2[#6][#6]1"
    ]
    
    has_linkage = False
    for pattern in linkage_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            has_linkage = True
            break

    # Check for characteristic hydroxylation patterns
    hydroxylation_patterns = [
        # A-ring patterns
        "Oc1cc(O)cc2c1",
        # B-ring patterns (including gallocatechin type)
        "c1c(O)c(O)ccc1",
        "c1c(O)c(O)c(O)cc1",
        # Gallate ester pattern
        "O=C(O)c1cc(O)c(O)c(O)c1"
    ]
    
    hydroxy_matches = 0
    for pattern in hydroxylation_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            hydroxy_matches += 1

    if hydroxy_matches < 2:
        return False, "Missing characteristic hydroxylation patterns"

    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = sum(1 for ring in ring_info.AtomRings() 
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))
    
    if aromatic_rings < 2:  # Need at least one A-ring and one B-ring
        return False, "Insufficient aromatic rings"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Lowered threshold
        return False, "Molecular weight too low for proanthocyanidin"

    # Final classification
    if total_flavan_matches >= 2 or (total_flavan_matches >= 1 and has_linkage):
        return True, "Contains flavan units with characteristic linkages and hydroxylation patterns"
    else:
        return False, "Does not meet minimum structural requirements for proanthocyanidin"