"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is defined as a molecule of high relative molecular mass,
    the structure of which essentially comprises the multiple repetition of units
    derived from molecules of low relative molecular mass.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Calculate molecular weight
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 1000:
            return False, f"Molecular weight is {mol_wt:.2f} Da, which is below typical macromolecule mass"

        # Define common linkage patterns in macromolecules
        peptide_bond = Chem.MolFromSmarts("C(=O)N")
        glycosidic_bond = Chem.MolFromSmarts("[CX4H]O[CX4H]")  # Simplified pattern for sugar linkages
        ester_bond = Chem.MolFromSmarts("C(=O)O[CX4H]")
        amide_bond = Chem.MolFromSmarts("C(=O)N[CX4H]")

        bond_patterns = [peptide_bond, glycosidic_bond, ester_bond, amide_bond]
        pattern_names = ['peptide bonds', 'glycosidic bonds', 'ester bonds', 'amide bonds']
        repeating_units_detected = []

        # Check for repeating linkage patterns
        for pattern, name in zip(bond_patterns, pattern_names):
            matches = mol.GetSubstructMatches(pattern)
            if len(matches) > 3:  # Threshold for considering repeating units
                repeating_units_detected.append(name)

        if not repeating_units_detected:
            return False, "No repeating units detected in the molecule"

        # If we reach here, the molecule is large and has repeating units
        reasons = f"Molecule is a macromolecule with molecular weight {mol_wt:.2f} Da and contains repeating units: {', '.join(repeating_units_detected)}"
        return True, reasons

    except Exception as e:
        return False, f"An error occurred during processing: {str(e)}"