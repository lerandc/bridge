import pathlib
import subprocess
from tqdm import tqdm
from ..io.utils import listfiles


def find_and_run_mult_edge(root_dir):
    # fps = [f for f in listfiles(root_dir) if "feff" in f and ".inp" in f]
    fps = [x for x in listfiles(root_dir)[0] if "inp" in x.name and "feff" in x.name]

    for f in tqdm(fps):
        # grab directory
        f_str = str(f)
        cur_dir = pathlib.Path(f_str.rsplit("/", 1)[0])
        cur_file = pathlib.Path(f_str)

        # grab edge and rename input to feff.inp
        cur_edge = f_str.rsplit("_", 1)[1][:-4]
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

        # mkdir for edge
        out_dir = pathlib.Path(cur_dir.joinpath(cur_edge))
        out_dir.mkdir(parents=True, exist_ok=True)

        # move xmu.dat to dir
        out_f = cur_dir.joinpath("xmu.dat")
        out_f.rename(out_dir.joinpath("xmu.dat"))

        # rename input back to original name
        cur_file_renamed.rename(cur_file)
