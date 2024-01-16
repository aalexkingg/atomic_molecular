from dataclasses import dataclass
import numpy as np
import scipy as sci
import pandas as pd

element_lookup = pd.read_csv("../base/Periodic Table of Elements.csv")


class Atom:
    def __init__(self, symbol):
        self.symbol: str

        self._lookup: pd.DataFrame = element_lookup.loc[element_lookup['Symbol'] == symbol]
        self.name: str = self._lookup["Element"]
        self.mass: float = self._lookup["AtomicMass"]     # amu
        self.charge: float = self._lookup["AtomicNumber"]   # C
        self.protons: int = self._lookup["NumberofProtons"]
        self.neutrons: int = self._lookup["NumberofNeutrons"]
        self.electrons: int = self._lookup["NumberofElectrons"]

    def __add__(self, other): ...
    def __sub__(self, other): ...


class Ion(Atom):
    def __int__(self, ion_charge):
        self.ion_charge: int = ion_charge
        super().charge += self.ion_charge
        self.Zcore = self.ion_charge + 1    # e



