from operator import itemgetter
from pymatgen.core import Structure, Molecule, Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.feff.outputs import LDos
from pymatgen.io.feff.inputs import  Potential, get_atom_map, get_absorbing_atom_symbol_index
from tabulate import tabulate


class MolPotential(Potential):
    """
    FEFF atomic potential for molecule or structure
    """

    def __init__(self, struct, absorbing_atom):
        """
        Args:
            struct (Structure): Structure or Molecule object.
            absorbing_atom (str/int): Absorbing atom symbol or site index
        """
        self.absorbing_atom, _ = get_absorbing_atom_symbol_index(absorbing_atom, struct)

        if struct.is_ordered:
            self.struct = struct
            self.pot_dict = get_atom_map(struct)
        else:
            raise ValueError("Structure with partial occupancies cannot be " "converted into atomic coordinates!")

        # if pulling into pymatgen later, probably should use base classes for IMolecule and IStructure for these checks
        if isinstance(self.struct, Structure):
            self.composition = self.struct.composition
        elif isinstance(self.struct, Molecule):
            # Adjust stoichiometries
            mol_copy = self.struct.copy()
            mol_copy.remove_sites([absorbing_atom])
            self.pot_dict = get_atom_map(mol_copy)
            self.composition = mol_copy.composition

    def __str__(self):
        """
        Returns a string representation of potential parameters to be used in
        the feff.inp file,
        determined from structure object.
                The lines are arranged as follows:
          ipot   Z   element   lmax1   lmax2   stoichiometry   spinph
        Returns:
            String representation of Atomic Coordinate Shells.
        """
        central_element = Element(self.absorbing_atom)

        if isinstance(self.struct, Structure):
            ipotrow = [[0, central_element.Z, central_element.symbol, -1, -1, 0.0001, 0]]
        elif isinstance(self.struct, Molecule):
            ipotrow = [[0, central_element.Z, central_element.symbol, -1, -1, 1, 0]]

        for el, amt in self.composition.items():
            ipot = self.pot_dict[el.symbol]
            ipotrow.append([ipot, el.Z, el.symbol, -1, -1, amt, 0])
        ipot_sorted = sorted(ipotrow, key=itemgetter(0))
        ipotrow = str(
            tabulate(
                ipot_sorted,
                headers=[
                    "*ipot",
                    "Z",
                    "tag",
                    "lmax1",
                    "lmax2",
                    "xnatph(stoichometry)",
                    "spinph",
                ],
            )
        )
        ipotlist = ipotrow.replace("--", "**")
        ipotlist = "".join(["POTENTIALS\n", ipotlist])

        return ipotlist