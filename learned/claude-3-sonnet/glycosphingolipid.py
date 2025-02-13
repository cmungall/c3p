"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: CHEBI:17670 glycosphingolipid
A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide. 
It is understood that the carbohydrate residue is attached by a glycosidic linkage to O-1 of the sphingoid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ceramide backbone (alkyl chain + amide + alkene)
    ceramide_pattern = Chem.MolFromSmarts("[N;X3][C;X3](=[O])[C;X4]([C;X4])(CCCC[C;X3]=C)")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"
    
    # Look for glycosidic linkage to sphingoid (O-C-O)
    glycosidic_pattern = Chem.MolFromSmarts("[O;X2][C;X4][O;X2]")
    if not AllChem.MolToSmarts(mol).count("OC") > 1:
        return False, "No glycosidic linkage found"
    
    # Look for carbohydrate residue (cyclitol or sugar)
    carb_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@@H]([C@H]([C@@H]1O)O)O)O")
    if not mol.HasSubstructMatch(carb_pattern):
        return False, "No carbohydrate residue found"
    
    # Check for multiple carbohydrate rings (glycosphingolipids often have branched oligosaccharides)
    rings = mol.GetRingInfo().AtomRings()
    carb_rings = [r for r in rings if all(mol.GetAtomWithIdx(i).GetSymbol() in ['C', 'O'] for i in r)]
    if len(carb_rings) < 2:
        return False, "Only one carbohydrate ring found, glycosphingolipids typically have branched oligosaccharides"
    
    # Check for common structural motifs of glycosphingolipids
    motifs = ['CC(=O)NC[C@H](O)/C=C', 'C[C@H](O)/C=C/C']
    mol_smarts = AllChem.MolToSmarts(mol)
    if not any(motif in mol_smarts for motif in motifs):
        return False, "No common structural motifs of glycosphingolipids found"
    
    # If all checks pass, classify as glycosphingolipid
    return True, "Contains a ceramide backbone with a carbohydrate residue attached via a glycosidic linkage"