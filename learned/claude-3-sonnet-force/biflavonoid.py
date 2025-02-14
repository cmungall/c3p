"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: CHEBI:33136 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer obtained by the oxidative coupling of at least two units
    of aryl-substituted benzopyran rings or its substituted derivatives, resulting in the two ring
    systems being joined together by a single atom or bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavonoid subunits
    flavonoid_pattern = Chem.MolFromSmarts("c1c(O)cc2c(c1O)C(=O)c3ccccc3O2")
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    
    # Check for at least 2 flavonoid subunits
    if len(flavonoid_matches) < 2:
        return False, "Molecule does not contain at least 2 flavonoid subunits"
    
    # Check that flavonoid subunits are joined by a single atom or bond
    joined_subunits = []
    for match1, match2 in itertools.combinations(flavonoid_matches, 2):
        match1_atoms = set(match1)
        match2_atoms = set(match2)
        
        # Check if subunits share a single atom
        shared_atoms = match1_atoms.intersection(match2_atoms)
        if len(shared_atoms) == 1:
            joined_subunits.append((match1, match2))
            continue
        
        # Check if subunits are joined by a bond
        for atom1 in match1_atoms:
            for atom2 in match2_atoms:
                if mol.GetBondBetweenAtoms(atom1, atom2):
                    joined_subunits.append((match1, match2))
                    break
    
    if len(joined_subunits) >= 1:
        return True, "Contains at least 2 flavonoid subunits joined by a single atom or bond"
    else:
        return False, "Flavonoid subunits not joined by a single atom or bond"