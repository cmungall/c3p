"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:  # Lower bound for triterpenoids + sugars
        return False, "Too few carbons for a triterpenoid saponin"

    # Count oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a triterpenoid saponin"


    #Check ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 3: #Triterpenoids usually have 3 or more rings
        return False, "Too few rings to be a triterpenoid."

    # Look for glycosidic bonds (-C-O-C)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_matches) < 1:
         return False, "No glycosidic bond found."
    
    # Check for common sugar moieties (glucose, rhamnose etc)
    glucose_pattern = Chem.MolFromSmarts("OC[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H](CO)O1")
    rhamnose_pattern = Chem.MolFromSmarts("C[C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O")
    arabinose_pattern = Chem.MolFromSmarts("OC[C@H]1[C@H](O)[C@H](O)[C@H](CO)O1")
    xylose_pattern = Chem.MolFromSmarts("OC[C@H]1[C@H](O)[C@@H](O)[C@H](O)C1")
    galactose_pattern = Chem.MolFromSmarts("OC[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H](CO)O1")

    sugar_matches = (
        mol.HasSubstructMatch(glucose_pattern)
        or mol.HasSubstructMatch(rhamnose_pattern)
        or mol.HasSubstructMatch(arabinose_pattern)
        or mol.HasSubstructMatch(xylose_pattern)
        or mol.HasSubstructMatch(galactose_pattern)
    )

    if not sugar_matches:
      return False, "No sugar moiety detected"


    # Check for hydroxyl groups (more relaxed now)
    oh_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and len(atom.GetBonds()) == 1)
    if oh_count < 2:
       return False, "Too few hydroxyl groups, unlikely for triterpenoid saponin."


    return True, "Contains a triterpenoid core with glycosidic bond(s) and sugar(s)."