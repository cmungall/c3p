"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:18374 phytosterol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    A phytosterol is a plant sterol similar to cholesterol, varying only in carbon side chains and/or presence or absence of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sterol core SMARTS pattern (four fused rings typical of sterols)
    sterol_core_smarts = 'C1CCC2(C1)C3CCC4(C(C3)C2)CCC=C4'  # Simplified sterol core
    sterol_core = Chem.MolFromSmarts(sterol_core_smarts)
    if sterol_core is None:
        return False, "Error in sterol core SMARTS pattern"

    # Check if molecule contains sterol core
    if not mol.HasSubstructMatch(sterol_core):
        return False, "Sterol core not found"

    # Generate Murcko scaffolds to compare the core structures
    mol_scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    mol_scaffold_smiles = Chem.MolToSmiles(mol_scaffold, isomericSmiles=True)

    # Cholesterol scaffold SMILES
    cholesterol_smiles = 'C[C@H](CCCC(C)C)C1CCC2C1(CC=C3C2CCC4=C3CC[C@@H](O)C4)C'  # Cholesterol
    cholesterol_mol = Chem.MolFromSmiles(cholesterol_smiles)
    cholesterol_scaffold = MurckoScaffold.GetScaffoldForMol(cholesterol_mol)
    cholesterol_scaffold_smiles = Chem.MolToSmiles(cholesterol_scaffold, isomericSmiles=True)

    # Check if the scaffolds are the same (allowing for stereochemistry differences)
    if mol_scaffold_smiles != cholesterol_scaffold_smiles:
        return False, "Molecule does not share the sterol scaffold with cholesterol"

    # Count the number of carbon atoms in side chains
    mol_atoms = mol.GetNumHeavyAtoms()
    scaffold_atoms = mol_scaffold.GetNumHeavyAtoms()
    side_chain_atoms = mol_atoms - scaffold_atoms

    # Phytosterols typically have longer side chains than cholesterol
    cholesterol_side_chain_atoms = cholesterol_mol.GetNumHeavyAtoms() - cholesterol_scaffold.GetNumHeavyAtoms()
    if side_chain_atoms < cholesterol_side_chain_atoms:
        return False, "Side chains are not consistent with phytosterols (too short)"

    # Check for differences only in side chains and/or double bonds
    # Remove side chains to compare the core
    mol_core = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts('*~[*]~[*]~[*]-[*]'))  # Remove side chains
    cholesterol_core = Chem.DeleteSubstructs(cholesterol_mol, Chem.MolFromSmarts('*~[*]~[*]~[*]-[*]'))

    # Check if cores are identical
    if not mol_core.HasSubstructMatch(cholesterol_core):
        return False, "Core structure differs beyond side chains and double bonds"

    return True, "Molecule contains sterol core and matches phytosterol criteria"