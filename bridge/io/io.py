import matplotlib.pyplot as plt
import pandas as pd
import pathlib
import pymatgen as mg
import re
import seaborn as sns
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin
from collections import OrderedDict
from reportlab.platypus import SimpleDocTemplate, Image
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch, cm
from reportlab.lib.pagesizes import letter, landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from .utils import listfiles

def plot_feff_edge(base_dir, target, save=False, target_dir="."):
    subfolders = [p for p in pathlib.Path(base_dir).iterdir() if "target_" + target in str(p)]
    if len(subfolders) == 0:
        return

    fig, axes = plt.subplots(2,1)
    avg_dict_x = OrderedDict()
    avg_dict_y = OrderedDict()
    num_sites = len(subfolders)
    for f in subfolders:
        edge_dirs = [p for p in f.iterdir() if (p.is_dir()) and not (p.name.startswith("."))]
        for i, e in enumerate(edge_dirs):
            edge = str(e).split("/")[-1]
            df = pd.read_csv(e.joinpath('xmu.dat'), skiprows=25,delimiter=r"\s+", \
                    names=['omega', 'e','k','mu','mu0','chi'])

            if edge in avg_dict_y.keys():
                avg_dict_y[edge] += df["mu"]/num_sites
            else:
                avg_dict_y[edge] = df["mu"]/num_sites
                avg_dict_x[edge+"_omega"] = df['omega']

            lbl = re.search("site_[0-9]+", str(f)).group(0).split('_')[1]
            axes[i].plot(df['omega'], df['mu'], color=(0.0,0.6,0.9),linewidth=.67)
            axes[i].text(0.9, 0.05, edge, transform=axes[i].transAxes)

    for i, e in enumerate(avg_dict_y.keys()):
        axes[i].plot(avg_dict_x[e+"_omega"], avg_dict_y[e], color=(0,0,0), linewidth=1)

    structure = (str(pathlib.Path(base_dir)).split("/")[-1]).rsplit("_",1)[0]
    axes[0].set_title(target + " edges for " + structure)
    axes[1].set_xlabel("omega")
    axes[0].set_ylabel("mu")
    axes[1].set_ylabel("mu")
    
    if save:
        p = pathlib.Path(target_dir)
        p.mkdir(parents=True, exist_ok=True)
        plt.savefig(p.joinpath(structure + "_xmu_" + target + "_edges.jpg"))
    else:
        plt.show()
    plt.close()


def plot_feff_absorber_dos(base_dir, target, save=False, target_dir="."):
    subfolders = [p for p in pathlib.Path(base_dir).iterdir() if "target_" + target in str(p)]
    if len(subfolders) == 0:
        return

    fig = plt.figure()
    ax = fig.gca()
    avg_dict = dict()
    num_sites = len(subfolders)
    for f in subfolders:
        df = pd.read_csv(f.joinpath('ldos00.dat'), skiprows=11,delimiter=r"\s+", 
                names=['e', 'sDOS','pDOS','dDOS','fDOS'])
        df['DOS'] = df.sDOS + df.pDOS + df.dDOS + df.fDOS
        if 'DOS' in avg_dict.keys():
            avg_dict['DOS'] += df['DOS'] / num_sites
        else:
            avg_dict['DOS'] = df['DOS'] / num_sites
            avg_dict['e'] = df['e']

        # lbl = re.search("site_[0-9]+", f.__str__()).group(0).split('_')[1]
        sns.lineplot(data=df, x='e',y='DOS',color=(0, 0.6, 0.9), linewidth=.67)

    # edge = f.__str__().split("_")[-1]
    ax.plot(avg_dict['e'], avg_dict['DOS'], color=(0,0,0), linewidth=1)
    structure = (pathlib.Path(base_dir).__str__().split("/")[-1]).rsplit("_",1)[0]
    ax.set_title("Feff Absorber DOS of " + target + " in " + structure)
    # plt.legend()

    if save:
        p = pathlib.Path(target_dir)
        p.mkdir(parents=True, exist_ok=True)
        plt.savefig(p.joinpath(structure+"_feff_dos_absorber_"+target+"_edges.jpg"))
    else:
        plt.show()

    plt.close()

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

def plot_feff_dos(base_dir, target, save=False, target_dir="."):

    subfolders = [p for p in pathlib.Path(base_dir).iterdir() if "target_" + target in str(p)]
    if len(subfolders) == 0:
        return

    try:
        pot_dict = get_pot_dict(subfolders[0].joinpath("feff.inp"))
    except FileNotFoundError:
        # search for an edge
        feff_inps = [f for f in subfolders[0].iterdir() \
                        if "inp" in f.name and "feff" in f.name]
        pot_dict = get_pot_dict(feff_inps[0])

    for elem in pot_dict.keys():
        fig = plt.figure()
        ax = fig.gca()
        avg_dict = dict()
        num_sites = len(subfolders)
        for f in subfolders:
            ldos_path = "ldos" + "%02i" % int(pot_dict[elem]) +".dat"
            df = pd.read_csv(f.joinpath(ldos_path), skiprows=11,delimiter=r"\s+", 
                    names=['e', 'sDOS','pDOS','dDOS','fDOS'])
            df['DOS'] = df.sDOS + df.pDOS + df.dDOS + df.fDOS

            if 'DOS' in avg_dict.keys():
                avg_dict['DOS'] += df['DOS'] / num_sites
            else:
                avg_dict['DOS'] = df['DOS'] / num_sites
                avg_dict['e'] = df['e']

            lbl = re.search("site_[0-9]+", str(f)).group(0).split('_')[1]
            sns.lineplot(data=df, x='e',y='DOS', color=(0, 0.6, 0.9), linewidth=.67)

        ax.plot(avg_dict['e'], avg_dict['DOS'], color=(0,0,0), linewidth=1)
        edge = str(f).split("_")[-1]
        structure = (str(pathlib.Path(base_dir)).split("/")[-1]).rsplit("_",1)[0]
        ax.set_title("Feff DOS for " + elem + " in " + structure)
        # plt.legend()

        if save:
            p = pathlib.Path(target_dir)
            p.mkdir(parents=True, exist_ok=True)
            plt.savefig(p.joinpath(structure + "_feff_dos_" + elem + "_" + target + "_edges.jpg"))
        else:
            plt.show()
        plt.close()

def plot_vasp_dos(directory, save=False, target_dir="."):
    dosrun = Vasprun(str(directory.joinpath('vasprun.xml')))

    species = set(dosrun.final_structure.species)

    indices = dict()
    for s in species:
        indices[s] = [i for i, x in enumerate(dosrun.final_structure.species) if x == s]
    
        fig = plt.figure()
        ax = fig.gca()
        num_sites = len(indices[s])
        for i, ind in enumerate(indices[s]):
            yup=dosrun.complete_dos.get_site_dos(dosrun.complete_dos.structure[ind]).densities[Spin.up]
            ydown = dosrun.complete_dos.get_site_dos(dosrun.complete_dos.structure[ind]).densities[Spin.down]
            y=yup+ydown
            if i == 0:
                y_avg = y/num_sites
            else:
                y_avg += y/num_sites
            x=dosrun.tdos.energies - dosrun.efermi
            ax.plot(x, y, label='{}'.format(str(i)), color=(0,0.6,0.9), linewidth=.67)
        
        ax.plot(x, y_avg, color=(0,0,0), linewidth=1)
        dir_tmp = str(pathlib.Path(directory)).rsplit("/",2)[-2:]
        sname = ""
        sname += dir_tmp[0] + "_" + dir_tmp[1]

        ax.set_title(s.name + " DOS for " +  sname + " calculated with VASP, PBE")
        ax.set_xlabel('e')
        ax.set_ylabel('DOS')
        # ax.set_xlim([-30, 30])
        # plt.legend()

        if save:
            p = pathlib.Path(target_dir)
            p.mkdir(parents=True, exist_ok=True)
            plt.savefig(p.joinpath(sname + "_vasp_dos_" + s.name + ".jpg"))
        else:
            plt.show()
        plt.close()


def make_report(fig_source_dir, outname):
    pdfmetrics.registerFont(TTFont("DejaVuSans", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"))

    all_files = sorted([x for x in listfiles(fig_source_dir)][0])

    prefix = re.split("_[0-9]+_[0-9]+_",str(all_files[0]).rsplit("/",1)[-1])[0]
    base_path = str(pathlib.Path(str(all_files[0]).rsplit("/",1)[0]).joinpath(prefix))

    structs = sorted(list(set([re.search("_[0-9]+_[0-9]+_",str(x).rsplit("/",1)[-1]).group(0)\
                                for x in all_files])))

    s_vasp_au = "vasp_dos_Au.jpg"
    s_vasp_p = "vasp_dos_P.jpg"
    s_feff_dos_au = "feff_dos_absorber_Au_edges.jpg"
    s_feff_dos_p = "feff_dos_absorber_P_edges.jpg"
    s_feff_xmu_au = "xmu_Au_edges.jpg"
    s_feff_xmu_p = "xmu_P_edges.jpg"

    cw = 10*inch
    ch = 5.625*inch
    c = canvas.Canvas(outname, pagesize=(cw, ch))
    c.setFont("DejaVuSans", size=20)
    for struct in structs:
        print("Adding page for " + struct[1:-1])
        #image position settings
        # 432W x 288H
        iw = 3*inch
        ih = iw*(288/432)
        wm = iw/10
        hm = ih/10
        x1 = wm
        x2 = iw + 2*wm
        x3 = 2*iw + 3*wm 
        y1 = hm
        y2 = ch - (ih+4*hm)

        #draw images
        # Au images
        try:
            c.drawImage(base_path+struct+s_vasp_au, x1, y2, iw, ih)
        except OSError:
            print("Failed VASP DOS Au for " + struct)

        try:
            c.drawImage(base_path+struct+s_feff_dos_au, x2, y2, iw, ih)
            c.drawImage(base_path+struct+s_feff_xmu_au, x3, y2, iw, ih)
        except OSError:
            print("Failed Feff outputs Au for " + struct)

        # P images
        if pathlib.Path(base_path+struct+s_feff_dos_p).exists():
            try:
                c.drawImage(base_path+struct+s_vasp_p, x1, y1, iw, ih)
            except OSError:
                print("Failed VASP DOS P for " + struct)

            try:
                c.drawImage(base_path+struct+s_feff_dos_p, x2, y1, iw, ih)
                c.drawImage(base_path+struct+s_feff_xmu_p, x3, y1, iw, ih)
            except OSError:
                print("Failed Feff outputs P for " + struct)

        #write titles and image subtitles
        c.setFontSize(size=20)
        c.drawCentredString(x=cw/2,y=ch-1*cm,text=struct[1:-1])

        c.showPage()
        
    c.save()