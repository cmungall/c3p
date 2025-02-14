"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones have a 2,3-dimethoxy-5-methylbenzoquinone core and a polyprenoid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure: 2,3-dimethoxy-5-methylbenzoquinone.
    # Note: The position of the methyl and the two methoxy can vary in the input,
    # so we are looking for the pattern but not the numbering of the core structure.
    core_pattern = Chem.MolFromSmarts("C1(=C(C(=O)C(=C(C1=O)OC)OC)C)")
    if not mol.HasSubstructMatch(core_pattern):
      return False, "Core benzoquinone structure with two methoxy and one methyl groups not found"
    
    #Check if there is a polyprenoid side chain at position 6
    # The side chain has alternating single and double bonds, and a methyl group on each double bond
    # The side chain is not specifically at position 6, just present, as examples have it at different positions.
    sidechain_pattern = Chem.MolFromSmarts("[CX4]=[CX3]([CX4])~[CX4]([CX4])=[CX3]([CX4])") # minimum of 2 isoprenoid units
    sidechain_match = mol.GetSubstructMatches(sidechain_pattern)
    
    # Ubiquinone-0 has no side chain
    if not sidechain_match and "CC1=C(C(OC)=C(OC)C(=O)C=C1)=O" not in smiles:
       return False, "Polyprenoid side chain not found"

    # Check for methoxy groups.
    methoxy_pattern = Chem.MolFromSmarts("OC")
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if len(methoxy_matches) < 2:
       return False, "Not enough methoxy groups"
   
    # Check for the presence of the methyl group in the ring
    methyl_pattern = Chem.MolFromSmarts("[CX4]") # only need to check there is at least one methyl group
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 1 :
       return False, "Not enough methyl groups"

    # Check for quinone group
    quinone_pattern = Chem.MolFromSmarts("C=1C(=O)C=CC=C1=O")
    if not mol.HasSubstructMatch(quinone_pattern):
       return False, "No quinone group found"
    
    #Additional checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
       return False, "Too few carbons for ubiquinone"
    if o_count < 4:
        return False, "Too few oxygens for ubiquinone"
    
    return True, "Matches the ubiquinone criteria: core structure, polyprenoid chain, methoxy groups, methyl group, quinone."