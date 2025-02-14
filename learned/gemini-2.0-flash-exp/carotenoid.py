"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are tetraterpenoids (C40) with a long conjugated polyene chain.
    Retinoids are excluded.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check for C40 Core Structure (approximate)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 38 or carbon_count > 42: # more relaxed count, typically 40
        if not mol.HasSubstructMatch(Chem.MolFromSmarts("C[C](C)(O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@H](CO)O1)")): #check for glycosylation
            return False, f"Incorrect carbon count: {carbon_count}. Must be close to 40 carbons for a basic carotenoid."

    # 3. Check for Long Conjugated Chain (at least 7 double bonds)
    # using SMARTS pattern to check for a conjugated system
    conjugated_chain_pattern = Chem.MolFromSmarts("[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]") #at least 7 double bonds
    if not mol.HasSubstructMatch(conjugated_chain_pattern):
        return False, "Does not contain a long conjugated chain of double bonds"

    # 4. Isoprenoid verification, look for isoprene patterns
    # This check is made optional
    # isoprene_pattern = Chem.MolFromSmarts("C[CH](C)[CH]=[CH]")
    # isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    # if len(isoprene_matches) < 4:  #Typically carotenoids have multiple isoprene units
    #     return False, "Does not have a typical isoprene pattern"

    # 5. Exclude Retinoids (beta-ionone ring, specific ring with double bond and methyls)
    retinoid_pattern = Chem.MolFromSmarts("CC1=C[CH](C)[CH2][CH2][C](C)=C1")
    if mol.HasSubstructMatch(retinoid_pattern):
         return False, "Contains a retinoid-like beta-ionone ring"

    # If all checks pass
    return True, "Meets criteria for a carotenoid: C40 core, long conjugated polyene, and no retinoid beta-ionone"