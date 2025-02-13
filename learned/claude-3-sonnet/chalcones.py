"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:35807 chalcone
A ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar, and its derivatives formed by substitution.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chalcone(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone contains two aromatic rings connected by an alpha,beta-unsaturated ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for chalcone backbone pattern: Ar-CH=CH-C(=O)-Ar
    chalcone_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)c2ccccc2")
    if not mol.HasSubstructMatch(chalcone_pattern):
        return False, "Missing chalcone backbone pattern: Ar-CH=CH-C(=O)-Ar"
    
    # Check for aromatic rings
    rings = mol.GetRingInfo().AtomRings()
    aromatic_rings = [ring for ring in rings if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"
    
    # Check for alpha,beta-unsaturated ketone
    ketone_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetDegree() == 1]
    for ketone_idx in ketone_atoms:
        # Look for C=C-C=O pattern
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, ketone_idx, 3)
        if len(env) == 4 and all(mol.GetBondBetweenAtoms(env[0], env[1]).GetBondType() == Chem.BondType.DOUBLE and
                                 mol.GetBondBetweenAtoms(env[1], env[2]).GetBondType() == Chem.BondType.SINGLE and
                                 mol.GetBondBetweenAtoms(env[2], env[3]).GetBondType() == Chem.BondType.DOUBLE):
            break  # Found alpha,beta-unsaturated ketone
    else:
        return False, "No alpha,beta-unsaturated ketone found"
    
    return True, "Contains two aromatic rings connected by an alpha,beta-unsaturated ketone (chalcone backbone)"