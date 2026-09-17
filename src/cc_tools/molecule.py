from dataclasses import dataclass

import numpy as np


@dataclass
class Molecule:
    symbols: list[str]
    coordinates: np.ndarray

    charge: int = 0
    multiplicity: int = 1

    name: str | None = None
    smiles: str | None = None

    def __post_init__(self):
        self.coordinates = np.asarray(self.coordinates, dtype=float)

        if self.coordinates.shape != (len(self.symbols), 3):
            raise ValueError(
                "coordinates must have shape (n_atoms, 3)"
            )