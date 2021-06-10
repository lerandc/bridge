import numpy as np
import pathlib
from pymatgen.core import Structure, Molecule, Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.feff.outputs import LDos
from pymatgen.io.feff.inputs import Atoms, Header, Potential, Tags

def poscar_to_xyz(fp, output="structure.xyz"):
    atoms_pos = Poscar.from_file(fp)
    mol = Molecule(atoms_pos.structure.species, atoms_pos.structure.cart_coords)
    mol.to("xyz", filename=output)

def poscar_to_mol(fp):
    atoms_pos = Poscar.from_file(fp)
    mol = Molecule(atoms_pos.structure.species, atoms_pos.structure.cart_coords)
    return mol

def create_feff_inp(atoms, target_species=None, target_dir="", **kwargs):
    # 3 options -> if atoms is molecule object; if atoms is poscar file; 


    if isinstance(atoms, Molecule):
        mol = atoms
    elif (atoms.split("/")[-1] == "POSCAR"):
        mol = poscar_to_mol(atoms)
    elif (atoms.split(".")[-1] == "xyz"):
        mol = Molecule.from_file(atoms)

    if target_species is None:
        target_species = mol.species[0]
    else:
        if not isinstance(target_species, Element):
            try:
                target_species = Element(target_species)
            except ValueError:
                print("Could not convert target species to element class.")

    # need a structure only to make header file so pymatgen can leader read the feff output
    lvecs = np.max(mol.cart_coords, axis=0)
    struct = mol.get_boxed_structure(lvecs[0], lvecs[1], lvecs[2])

    tags = {"CONTROL": "1 1 1 1 1 1",
            "RPATH": -1,
            "SCF": "7 0 30 0.2 3",
            "FMS": "9 0",
            "LDOS": "-30.0 30.0 0.1",
            "EDGE": "K",
            "COREHOLE": "RPA",
            "XANES":""
            }

    for key in kwargs:
        tags[key] = kwargs[key]

    fheader = Header(struct)
    ftags = Tags(tags)

    # TODO: reduce over unique site indices-- is there a pymatgen routine for molecule site equivalence?
    # then make subdirectory in target_dir for each site, named "target_Z_site_N_edge_K"
    ind = [i for i, j in enumerate(mol.species) if j is target_species]

    for i, j in enumerate(ind):
        fpot = Potential(mol, absorbing_atom=j)
        fatoms = Atoms(mol, absorbing_atom=j, radius=50)

        out_str = fheader.__str__() + "\n\n" \
                    + ftags.__str__() + "\n\n" \
                    + fpot.__str__() + "\n\n" \
                    + fatoms.__str__()

        out_dir = target_dir + "target_" + target_species.name \
                    + "_site_" + str(i) + "_edge_" + tags["EDGE"] + "/"
                
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
        with open(out_dir+"feff.inp", "w+") as f:
            f.write(out_str)
            f.close()
