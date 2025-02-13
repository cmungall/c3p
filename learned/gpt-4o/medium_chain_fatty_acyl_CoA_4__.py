"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA moiety
    coa_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H]([C@@H]1OP([O-])([O-])=O)OCN2CN=C(N)N=C2")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Check for medium-chain fatty acyl group, length 6 to 12 carbons
    # Note: Define pattern for medium chain fatty acyl part
    fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No fatty acyl chain connected via thioester bond found"

    # Extract the aliphatic chain preceding the thioester bond
    chain_carbon_atoms = []

    for atom in mol.GetAtoms():
        # Only consider aliphatic carbon atoms
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 4:
            chain_carbon_atoms.append(atom)
    
    # Check the length of the chain is within 6 to 12 carbons
    if len(chain_carbon_atoms) < 6 or len(chain_carbon_atoms) > 12:
        return False, f"Carbon chain length {len(chain_carbon_atoms)} is not medium-chain (6-12 carbons)"

    # Check deprotonation (4- charge) of phosphate and diphosphate groups
    phosphate_matches = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1]
    if len(phosphate_matches) < 4:
        return False, f"Insufficient deprotonated oxygens detected, found {len(phosphate_matches)}"

    return True, "Molecule is a medium-chain fatty acyl-CoA(4-)"