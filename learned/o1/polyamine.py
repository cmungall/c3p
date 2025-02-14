"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:18154 polyamine
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    Amino groups include primary, secondary, and tertiary amines (both aliphatic and aromatic),
    including their protonated forms, but excluding nitrogen atoms in amides, nitro groups,
    nitriles, azo compounds, and other non-amino functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    amino_nitrogens = set()

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:  # Not a nitrogen atom
            continue

        # Initialize exclusion flags
        is_amide = False
        is_nitro = False
        is_nitrile = False
        is_azo = False
        is_quaternary = False
        is_in_ring = atom.IsInRing()
        is_aniline = False

        # Exclude quaternary ammonium (N+ with 4 bonds)
        if atom.GetFormalCharge() > 0 and atom.GetTotalDegree() >= 4:
            is_quaternary = True

        if is_quaternary:
            continue

        # Check if nitrogen is in amide (N-C(=O))
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 6 and bond.GetBondType() == rdchem.BondType.SINGLE:
                # Check if carbon has a double bond to oxygen
                for nbr_bond in nbr.GetBonds():
                    nbr_atom = nbr_bond.GetOtherAtom(nbr)
                    if (nbr_atom.GetAtomicNum() == 8 and
                        nbr_bond.GetBondType() == rdchem.BondType.DOUBLE):
                        is_amide = True
                        break
                if is_amide:
                    break

        if is_amide:
            continue

        # Check if nitrogen is in nitro group (N(=O)-O)
        oxy_double = False
        oxy_single = False
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    oxy_double = True
                elif bond.GetBondType() == rdchem.BondType.SINGLE:
                    oxy_single = True
        if oxy_double and oxy_single:
            is_nitro = True

        if is_nitro:
            continue

        # Check if nitrogen is in nitrile (N#C)
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.TRIPLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 6:
                    is_nitrile = True
                    break

        if is_nitrile:
            continue

        # Check if nitrogen is in azo group (N=N)
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 7:
                    is_azo = True
                    break

        if is_azo:
            continue

        # Exclude nitrogen atoms in rings unless they are aniline-type amines
        if is_in_ring:
            # Check if nitrogen is attached to an aromatic carbon outside the ring (aniline)
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 6:
                    if nbr.GetIsAromatic() and not nbr.IsInRing():
                        is_aniline = True
                        break
            if not is_aniline:
                continue  # Exclude nitrogen atoms in rings not aniline-type

        # Include protonated amines (e.g., [NH3+])
        if (atom.GetHybridization() == rdchem.HybridizationType.SP3 or atom.GetIsAromatic()) and not is_quaternary:
            amino_nitrogens.add(atom.GetIdx())

    amino_count = len(amino_nitrogens)

    if amino_count >= 2:
        return True, f"Molecule contains {amino_count} amino groups"
    else:
        return False, f"Molecule contains {amino_count} amino group(s), need at least 2"