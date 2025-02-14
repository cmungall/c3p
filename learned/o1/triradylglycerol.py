"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol is a glycerol compound where each of the three hydroxyl groups 
    is substituted with either an acyl, alkyl, or alk-1-enyl group at positions 
    sn-1, sn-2, and sn-3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the glycerol backbone where each carbon is connected to an oxygen
    glycerol_pattern = Chem.MolFromSmarts("[C]-[C]-[C]")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Iterate through the possible glycerol backbones
    for match in matches:
        c1_idx, c2_idx, c3_idx = match
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Check that each carbon has an oxygen connected via single bond
        c_atoms = [c1, c2, c3]
        valid = True
        reasons = []
        for idx, c_atom in enumerate(c_atoms):
            oxygen_found = False
            for neighbor in c_atom.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), neighbor.GetIdx())
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Found an oxygen connected via single bond
                    oxygen = neighbor
                    # Check if the oxygen is connected to a substituent (acyl, alkyl, alk-1-enyl)
                    substituents = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetIdx() != c_atom.GetIdx()]
                    if not substituents:
                        valid = False
                        reasons.append(f"Oxygen at position {idx+1} has no substituent")
                        break
                    else:
                        substituent_atom = substituents[0]
                        # Check for acyl group (O=C-C), alkyl group (C), or alk-1-enyl group (C=C-C)
                        acyl_pattern = Chem.MolFromSmarts("C(=O)[#6]")
                        alkyl_pattern = Chem.MolFromSmarts("[#6]")
                        alk1enyl_pattern = Chem.MolFromSmarts("C=C[#6]")
                        # Create the fragment starting from substituent_atom
                        env = Chem.FindAtomEnvironmentOfRadiusN(mol, 2, substituent_atom.GetIdx())
                        amap = {}
                        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
                        if submol.HasSubstructMatch(acyl_pattern):
                            reasons.append(f"Position {idx+1}: acyl substituent found")
                            oxygen_found = True
                        elif submol.HasSubstructMatch(alk1enyl_pattern):
                            reasons.append(f"Position {idx+1}: alk-1-enyl substituent found")
                            oxygen_found = True
                        elif submol.HasSubstructMatch(alkyl_pattern):
                            reasons.append(f"Position {idx+1}: alkyl substituent found")
                            oxygen_found = True
                        else:
                            valid = False
                            reasons.append(f"Position {idx+1}: substituent is not acyl, alkyl, or alk-1-enyl")
                            break
                elif neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    # Possibly a carbonyl oxygen, skip
                    continue
            if not oxygen_found:
                valid = False
                reasons.append(f"Carbon at position {idx+1} is not connected to an appropriate oxygen")
                break
        if valid:
            return True, "Molecule is a triradylglycerol"
        else:
            continue  # Try next possible glycerol backbone

    return False, "; ".join(reasons) if reasons else "No valid triradylglycerol structure found"