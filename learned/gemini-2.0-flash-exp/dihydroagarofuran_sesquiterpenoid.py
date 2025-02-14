"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid has a specific bicyclic core with an ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
                         Reason for the classification
    """

    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Core Skeleton Matching using SMARTS
    # This SMARTS matches the dihydroagarofuran core, with * as variable carbon positions.
    # The stereochemistry is somewhat included, but not fully. Note that the 
    # methyls are not defined here, they are added later.
    core_pattern = Chem.MolFromSmarts('[C@]([H])1([C@]([H])2[C@]([H])([C]([H])([H])[O]1)[C]([H])([H])[C]([H])([H])[C]([H])([C]2([H])[H])[H])')
    if not mol.HasSubstructMatch(core_pattern):
         return False, "No dihydroagarofuran core skeleton found."

    # 3. Check for at least 1 methyl group on the core
    methyl_pattern = Chem.MolFromSmarts('[C]([H])([H])([H])')
    core_match = mol.GetSubstructMatch(core_pattern)
    if core_match:
        methyl_count = 0
        for atom_idx in core_match:
          atom = mol.GetAtomWithIdx(atom_idx)
          for neighbor in atom.GetNeighbors():
            if mol.HasSubstructMatch(methyl_pattern, neighbor.GetIdx()):
                 methyl_count+=1
        if methyl_count <1 :
            return False, "No methyl group directly bonded to core carbon"
    else:
        return False, "No core pattern found, should not reach here"
    
    # 4. Additional checks: At least 2 ester groups and 4 oxygens
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
         return False, f"Less than 2 ester groups. Found {len(ester_matches)}"
    
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 4:
        return False, f"Too few oxygens. Found {oxygen_count}"
    
    
    # 5. If all checks pass, it's classified as a dihydroagarofuran sesquiterpenoid.
    return True, "Contains a dihydroagarofuran core and at least one methyl group on the core, at least two ester groups and four oxygens"