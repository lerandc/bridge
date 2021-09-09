import pathlib
import subprocess
from functools import reduce
from tqdm import tqdm
from ..io.utils import listfiles


def find_and_run_mult_edge(root_dir, rerun=False):
    # fps = [f for f in listfiles(root_dir) if "feff" in f and ".inp" in f]
    fps = [x for x in listfiles(root_dir) if "inp" in x.name and "feff" in x.name]

    for f in tqdm(fps):
        # grab directory
        f_str = str(f)
        cur_dir = pathlib.Path(f_str.rsplit("/", 1)[0])
        cur_file = pathlib.Path(f_str)

        # grab edge and rename input to feff.inp
        # expect filenames of feff_****.inp where **** indicates edge/sim type for structure
        cur_edge = f_str.rsplit("/", 1)[1].split("_", 1)[1][:-4]
        
        with open(cur_file, "r") as f:
                lines = f.readlines()
                rfunc = lambda x, y: x or y
                elnes_check = reduce(rfunc, ["ELNES" in x for x in lines])
                exelfs_check = reduce(rfunc, ["EXELFS" in x for x in lines])

                is_eels = elnes_check or exelfs_check

        if is_eels:
            out_f_str = "eels.dat"
        else:
            out_f_str = "xmu.dat"

        out_dir = pathlib.Path(cur_dir.joinpath(cur_edge))

        output_exists = pathlib.Path(out_dir.joinpath(out_f_str)).is_file()
        if rerun or (not output_exists):
            cur_file.rename(cur_dir.joinpath("feff.inp"))
            cur_file_renamed = pathlib.Path(cur_dir.joinpath("feff.inp"))

            # check if pot bin exists in dir
            pot_f = cur_dir.joinpath("pot.bin")

            if pot_f.is_file():
                # run feff_new_edge
                process = "feff_new_edge"
            else:
                # run feff
                process = "feff"

            f_stdout = open(cur_dir.joinpath("stdout"), "w")
            f_stderr = open(cur_dir.joinpath("stderr"), "w")
            subprocess.call(process, cwd=cur_dir, stdout=f_stdout, stderr=f_stderr)


            # check for existence of output
            out_f = cur_dir.joinpath(out_f_str)
            if out_f.is_file():
                # mkdir for edge
                out_dir.mkdir(parents=True, exist_ok=True)

                # move xmu.dat to dir
                out_f.rename(out_dir.joinpath(out_f_str))
            else:
                print("Simulation " + str(f) + " failed. Renaming file and moving on.")

            # rename input back to original name
            cur_file_renamed.rename(cur_file)
