"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a mineral, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of charged atoms
    has_charge = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())
    if not has_charge:
        return False, "No charged atoms, unlikely to be mineral"

    # Check for metal presence
    metal_pattern = Chem.MolFromSmarts("[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Y,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Tc,Re,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Sn,Pb,Sb,Bi,Po,Ge,As,Sb,Se,Te,B,Si]")
    if not mol.HasSubstructMatch(metal_pattern):
        return False, "No metals detected, unlikely to be mineral"

    # Check for common anions (oxide, hydroxide, carbonate, sulfate, phosphate, halides)
    anion_patterns = [
        Chem.MolFromSmarts("[O-]"),      # oxide
        Chem.MolFromSmarts("[OH-]"),    # hydroxide
        Chem.MolFromSmarts("C(=O)([O-])[O-]"),   # carbonate
        Chem.MolFromSmarts("S(=O)(=O)([O-])[O-]"), # sulfate
        Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]"), # phosphate
        Chem.MolFromSmarts("[F-]"),      # fluoride
        Chem.MolFromSmarts("[Cl-]"),    # chloride
        Chem.MolFromSmarts("[Br-]"),    # bromide
        Chem.MolFromSmarts("[I-]"),      # iodide

    ]
    has_anion = any(mol.HasSubstructMatch(pattern) for pattern in anion_patterns)
    if not has_anion:
        return False, "No common mineral anions detected"

    #check charge balance, it's OK if it's an integer
    charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if abs(charge) > 2:
        return False, "Net charge is not likely for a mineral"

    # Check for water molecules if there are any separate water molecules
    water_pattern = Chem.MolFromSmarts("[O]")
    water_count = len(mol.GetSubstructMatches(water_pattern))
    if water_count > 0 :
      # If there are any 'O' atoms, and there are also charged species it could be a hydrated mineral
      if not has_charge:
          return False, "Contains neutral water but no charged species."


    # Check for long carbon chains
    carbon_chain_pattern = Chem.MolFromSmarts("[C]([C])[C]") #detect 3+ carbon chain
    if mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Contains carbon chains, unlikely to be mineral"


    return True, "Likely a mineral (contains metal, common anions and is charge balanced)"