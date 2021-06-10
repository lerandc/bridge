import numpy as np
from pymatgen.core import Structure, Molecule, Lattice
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.feff.outputs import LDos
from pymatgen.io.feff.inputs import Atoms, Header, Potential, Tags
from atomate.feff.workflows.core import get_unique_site_indices

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

    # need a structure only to make header file so pymatgen can leader read the feff output
    lvecs = np.max(mol.cart_coords, axis=0)
    struct = mol.get_boxed_structure(lvecs[0], lvecs[1], lvecs[2])

    tags = {"CONTROL": "1 1 1 1 1 1",
            "RPATH": -1,
            "SCF": "7 0 30 0.2 3",
            "FMS": "9 0",
            "LDOS": "-30.0 30.0 0.1",
            "RECIPROCAL": "",
            "EDGE": "K",
            "COREHOLE": "RPA",
            }

    for key in kwargs:
        tags[key] = kwargs[key]

    fheader = Header(struct)
    ftags = Tags(tags)

    # TODO: reduce over unique site indices
    # then make subdirectory in target_dir for each site 

    fpot = Potential(mol, absorbing_atom=0)
    fatoms = Atoms(mol, absorbing_atom=0, radius=50)

    out_str = fheader.__str__() + "\n\n" \
                + ftags.__str__() + "\n\n" \
                + fpot.__str__() + "\n\n" \
                + fatoms.__str__()
            

    with open(target_dir+"feff.inp", "w+") as f:
        f.write(out_str)
        f.close()
