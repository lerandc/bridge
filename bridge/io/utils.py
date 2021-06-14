import numpy as np
import pathlib
import os
from pymatgen.core import Structure, Molecule, Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.feff.outputs import LDos
from pymatgen.io.feff.inputs import Header, Tags
from .feff import MolPotential, MolAtoms

def poscar_to_xyz(fp, output="structure.xyz"):
    atoms_pos = Poscar.from_file(fp)
    mol = Molecule(atoms_pos.structure.species, atoms_pos.structure.cart_coords)
    mol.to("xyz", filename=output)

def poscar_to_mol(fp):
    atoms_pos = Poscar.from_file(fp)
    mol = Molecule(atoms_pos.structure.species, atoms_pos.structure.cart_coords)
    return mol

def create_feff_input(atoms, target_species=None, target_dir="", **kwargs):
    # 3 options -> if atoms is molecule object; if atoms is poscar file; 

    if isinstance(atoms, Molecule):
        mol = atoms
    elif (atoms.split("/")[-1] == "POSCAR"):
        mol = poscar_to_mol(atoms)
    elif (atoms.split("/")[-1] == "CONTCAR"):
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
    ind = [i for i, j in enumerate(mol.species) if j is target_species]

    root_path = pathlib.Path(target_dir)
    for i, j in enumerate(ind):
        fpot = MolPotential(mol, absorbing_atom=j)
        fatoms = MolAtoms(mol, absorbing_atom=j, radius=50)

        # add overlap correction for hydrogen as per feff guidebook pg 96/97
        if "H" in fpot.pot_dict:
            folp_string = "FOLP " + str(fpot.pot_dict["H"]) + " 0.8\n"
        else:
            folp_string = ""

        out_str = fheader.__str__() + "\n\n" \
                    + ftags.__str__() + "\n" \
                    + folp_string + "\n" \
                    + fpot.__str__() + "\n\n" \
                    + fatoms.__str__()

        out_dir = root_path.joinpath("target_" + target_species.name \
                    + "_site_" + str(i) + "_edge_" + tags["EDGE"])
                
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
        with open(out_dir.joinpath("feff.inp"), "w+") as f:
            f.write(out_str)
            f.close()

def listfiles(folder):
    # generator for files in subdirectory
    # https://stackoverflow.com/questions/12420779/simplest-way-to-get-the-equivalent-of-find-in-python
    for root, folders, files in os.walk(folder):
        for filename in folders + files:
            yield os.path.join(root, filename)

def find_vasp_outputs_make_inputs(source_dir, target_root_dir, target_species, file="CONTCAR", **kwargs):
    fps = [f for f in listfiles(source_dir) if f.split("/")[-1]==file]

    root_path = pathlib.Path(target_root_dir)
    for f in fps:
        full_path = f.split("/")
        sub_dir = ""
        for folder in full_path[:-1]:
            sub_dir += folder + "_"
        
        sub_dir += "feff"

        target_path = root_path.joinpath(sub_dir)

        for t in target_species:
            create_feff_input(f, t, target_path, **kwargs)