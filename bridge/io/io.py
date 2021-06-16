import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import pymatgen as mg
import re
import seaborn as sns
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin

def plot_feff_edge(base_dir, target):
    subfolders = [p for p in pathlib.Path(base_dir).iterdir() if "target_" + target in p.__str__()]
    for f in subfolders:
        df = pd.read_csv(f.joinpath('xmu.dat'), skiprows=25,delimiter=r"\s+", 
                names=['omega', 'e','k','mu','mu0','chi'])

        lbl = re.search("site_[0-9]+_", f.__str__()).group(0).split('_')[1]
        sns.lineplot(data=df, x='omega',y='mu',label=target+" #"+lbl)

    edge = f.__str__().split("_")[-1]
    structure = (pathlib.Path(base_dir).__str__().split("/")[-1]).rsplit("_",1)[0]
    plt.title(target + " " + edge + " edge for " + structure)
    plt.legend()
    plt.show()


def plot_feff_absorber_dos(base_dir, target):
    subfolders = [p for p in pathlib.Path(base_dir).iterdir() if "target_" + target in p.__str__()]
    for f in subfolders:
        df = pd.read_csv(f.joinpath('ldos00.dat'), skiprows=11,delimiter=r"\s+", 
                names=['e', 'sDOS','pDOS','dDOS','fDOS'])
        df['DOS'] = df.sDOS + df.pDOS + df.dDOS + df.fDOS
        lbl = re.search("site_[0-9]+_", f.__str__()).group(0).split('_')[1]
        sns.lineplot(data=df, x='e',y='DOS',label=target+" #"+lbl)

    edge = f.__str__().split("_")[-1]
    structure = (pathlib.Path(base_dir).__str__().split("/")[-1]).rsplit("_",1)[0]
    plt.title("Absorber DOS of " + target + " " + edge + " edge for " + structure)
    plt.legend()
    plt.show()

def get_pot_dict(fp):
    """ fp is path to feff.inp"""
    with open(fp, "r") as f:
        for i,line in enumerate(f):
            if line == "POTENTIALS\n":
                break
        for j, line in enumerate(f):
            if line == "ATOMS\n":
                break
    f = open(fp,"r")
    ipot_lines = f.readlines()[i:i+j]
    ipot_lines = [i.split() for i in ipot_lines[3:]]
    pot_dict = dict()
    for pot in ipot_lines:
        pot_dict[pot[2]] = pot[0]

    return pot_dict

def plot_feff_dos(base_dir, target, elem):

    subfolders = [p for p in pathlib.Path(base_dir).iterdir() if "target_" + target in p.__str__()]

    pot_dict = get_pot_dict(subfolders[0].joinpath("feff.inp"))
    for f in subfolders:
        ldos_path = "ldos" + "%02i" % int(pot_dict[elem]) +".dat"
        df = pd.read_csv(f.joinpath(ldos_path), skiprows=11,delimiter=r"\s+", 
                names=['e', 'sDOS','pDOS','dDOS','fDOS'])
        df['DOS'] = df.sDOS + df.pDOS + df.dDOS + df.fDOS
        lbl = re.search("site_[0-9]+_", f.__str__()).group(0).split('_')[1]
        sns.lineplot(data=df, x='e',y='DOS',label=target+" #"+lbl)

    edge = f.__str__().split("_")[-1]
    structure = (pathlib.Path(base_dir).__str__().split("/")[-1]).rsplit("_",1)[0]
    plt.title(elem + " DOS of " + target + " " + edge + " edge for " + structure)
    plt.legend()
    plt.show()

def plot_vasp_dos(directory):
    dosrun = Vasprun(directory+'vasprun.xml')

    species = set(dosrun.final_structure.species)

    indices = dict()
    for s in species:
        indices[s] = [i for i, x in enumerate(dosrun.final_structure.species) if x == s]
    
        for i in indices[s]:
            yup=dosrun.complete_dos.get_site_dos(dosrun.complete_dos.structure[i]).densities[Spin.up]
            ydown = dosrun.complete_dos.get_site_dos(dosrun.complete_dos.structure[i]).densities[Spin.down]
            y=yup+ydown
            x=dosrun.tdos.energies - dosrun.efermi
            plt.plot(x, y, label='{}'.format(str(i)))
        
        sname = ""
        for i in directory.rsplit("/",2)[-2:]:
            sname+=i

        plt.title(s.name + " DOS for " +  sname + " calculated with VASP, PBE")
        plt.xlabel('e')
        plt.ylabel('DOS')
        plt.legend()
        plt.show()