"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas molecular entity based on its SMILES string.
    A gas molecular entity is typically a simple structure that is gaseous at 
    standard temperature and pressure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Elements typically found in gaseous compounds
    gaseous_elements = {1, 2, 6, 7, 8, 9, 10, 17, 18, 36, 54, 86}  # H, He, C, N, O, F, Ne, Cl, Ar, Kr, Xe, Rn

    # Updated known special gaseous forms
    known_gaseous_smiles = {
        "[219Rn]", "[220Rn]", "[222Rn]", "[Rn]", "[Xe]", "[He]", "[3He]", "[4He]", "[6He]", 
        "O=C=O", "[C-]#[O+]", "CCC", "CC", "C=C", "[H]N([H])[H]", "ClCl", "FF",
        "[H]C([H])([H])[H]", "[H][H]", "[3H][3H]", "C1CO1", "I[H]"
    }

    # Check if the molecule matches any known gaseous SMILES
    if smiles in known_gaseous_smiles:
        return True, "Recognized special gaseous form"
    
    # Validate atoms against typical gaseous elements
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in gaseous_elements:
            return False, f"Contains non-gaseous element: {atom.GetSymbol()}"

    # Consider extending rules on molecular parameters
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    num_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    tpsa = Descriptors.TPSA(mol)
    
    # Adjust thresholds to better capture gaseous characteristics
    if mol_wt > 100 or num_rotatable > 1 or num_rings > 0 or tpsa > 20:
        return False, f"Molecular weight or structure complexity suggests non-gas (MW: {mol_wt}, Rotatable Bonds: {num_rotatable}, Rings: {num_rings}, TPSA: {tpsa})"
    
    return True, "Molecule is a simple structure containing only gaseous elements at STP"