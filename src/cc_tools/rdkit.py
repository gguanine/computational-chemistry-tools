from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from .molecule import Molecule


def from_rdkit(
    mol: Chem.Mol,
    *,
    name: str | None = None,
    smiles: str | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
) -> Molecule:
    conformer = mol.GetConformer()

    symbols = [
        atom.GetSymbol()
        for atom in mol.GetAtoms()
    ]

    coordinates = [
        list(conformer.GetAtomPosition(i))
        for i in range(mol.GetNumAtoms())
    ]

    if charge is None:
        charge = Chem.GetFormalCharge(mol)

    if multiplicity is None:
        multiplicity = Descriptors.NumRadicalElectrons(mol) + 1

    if name is None:
        name = Chem.inchi.MolToInchiKey(mol)

    return Molecule(
        symbols=symbols,
        coordinates=coordinates,
        charge=charge,
        multiplicity=multiplicity,
        name=name,
        smiles=smiles,
    )


def from_smiles(
    smiles: str,
    *,
    name: str | None = None,
    charge: int | None = None,
    multiplicity: int | None = None,
    random_seed: int = 6,
) -> Molecule:
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    result = AllChem.EmbedMolecule(
        mol,
        randomSeed=random_seed,
    )

    if result != 0:
        raise RuntimeError("RDKit 3D embedding failed")

    return from_rdkit(
        mol,
        name=name,
        smiles=smiles,
        charge=charge,
        multiplicity=multiplicity,
    )