"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    An icosanoid is a signaling molecule arising from oxidation of C20 essential fatty acids
    such as EPA, AA, and DGLA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count number of carbon atoms
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(c_atoms)
    if num_carbons < 18 or num_carbons > 22:
        return False, f"Number of carbon atoms ({num_carbons}) not in range for icosanoids (18-22)"
    
    # Check for oxygen-containing functional groups
    o_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(o_atoms) < 2:
        return False, "Less than 2 oxygen atoms found, unlikely to be an oxidized fatty acid"
    
    # Check for carboxylic acid group (-COOH)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found, common in icosanoids"
    
    # Check for number of double bonds
    num_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            num_double_bonds += 1
    if num_double_bonds < 2:
        return False, f"Only {num_double_bonds} double bonds found, less than expected for icosanoids"
    
    # Optional: Check for specific functional groups (epoxy, keto)
    keto_pattern = Chem.MolFromSmarts("C(=O)")
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    has_keto = mol.HasSubstructMatch(keto_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    if not (has_keto or has_epoxy):
        return False, "No keto or epoxy groups found, common in icosanoids"
    
    # Estimate if molecule could be derived from C20 fatty acids
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 280 or mol_wt > 380:
        return False, f"Molecular weight ({mol_wt:.2f}) outside typical range for icosanoids"
    
    # Lipinski's Rule of Five compliance (optional, icosanoids are bioactive lipids)
    if Lipinski.NumHDonors(mol) > 5 or Lipinski.NumHAcceptors(mol) > 10:
        return False, "Does not comply with hydrogen bond donor/acceptor count for typical icosanoids"
    
    return True, "Molecule matches key features of icosanoids: C20 backbone, oxidation products, multiple double bonds"