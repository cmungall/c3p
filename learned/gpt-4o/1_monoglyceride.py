"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride contains a glycerol backbone with one fatty acid chain
    attached via an ester bond at the 1-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Glycerol backbone with ester linkage at 1-position
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    ester_pattern = Chem.MolFromSmarts("O(C=O)")
    
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester linkage, found {len(ester_matches)}"
    
    # Check position of ester linkage
    ester_sites = [bond.GetBeginAtom().GetIdx() for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.ESTER]
    if not ester_sites or ester_sites[0] != 0:
        return False, "Ester linkage not at 1-position"

    # Verify the presence of a single acyl chain attached at the ester group
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)C")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_chain_matches) != 1:
        return False, "Expected 1 acyl chain, found {}".format(len(acyl_chain_matches))
    
    # Consideration of chiral center
    n_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=False))
    if n_chiral_centers != 1:
        has_chiral = "chiral" if n_chiral_centers >= 1 else "achiral"
        return False, f"Molecule should ideally have 1 chiral center, found {n_chiral_centers} ({has_chiral})"
    
    return True, "Valid structure for a 1-monoglyceride"