#!/usr/bin/env python

import os
import re
import sys
from pathlib import Path
from shutil import move
from subprocess import CompletedProcess
from typing import Optional

from tadashi import TrEnum
from tadashi.apps import App
from tadashi.translators import Pet, Polly, Translator

ml4tadashi = os.path.dirname(__file__)
ml4tadashi = os.path.dirname(ml4tadashi)
ml4tadashi = os.path.dirname(ml4tadashi)
ml4tadashi = os.path.dirname(ml4tadashi)
ml4tadashi = os.path.dirname(ml4tadashi)
sys.path.append(ml4tadashi)
ml4tadashi = os.path.join(ml4tadashi, "ML4TADASHI")
sys.path.append(ml4tadashi)

import ML4TADASHI


class CloverLeaf(App):
    def __init__(
        self,
        source="kernels/advec_mom_kernel_c.c",
        translator: Optional[Translator] = None,
        compiler_options: list = None,
        ephemeral: bool = False,
        populate_scops: bool = True,
    ):
        super().__init__(
            source=source,
            translator=translator,
            compiler_options=compiler_options,
            ephemeral=ephemeral,
            populate_scops=populate_scops,
        )

    def codegen_init_args(self):
        return {}

    def app_required_options(self) -> list[str]:
        return []

    def compile_cmd(self, suffix: str) -> list[str]:
        src = self.source.with_suffix(".o")
        if src.exists() and isinstance(self.translator, Polly):
            dst = src.parents[1] / src.name
            move(src=src, dst=dst)
            src.touch()
        cmd = [
            "make",
            "-j",
            f"SOURCE={self.source.with_suffix('').name}",
            # f"CSOURCE={src.with_suffix('').name}",
            "COMPILER=GNU",
            "C_MPI_COMPILER=mpicc",
            "MPI_COMPILER=mpif90",
        ]
        return cmd

    @property
    def output_binary(self) -> Path:
        outbin = self.source.parents[1] / self.source.with_suffix(".x").name
        outbin = Path(f"./{outbin}")
        outbin.is_absolute()
        return outbin

    def extract_runtime(self, proc: CompletedProcess) -> float:
        stdout = proc.stderr.decode()
        lines = stdout.split("\n")
        pattern = re.compile(r".*Wall clock +(\d+\.\d+).*")
        rv = -1.0
        for line in lines:
            match = pattern.match(line)
            if match:
                rv = float(match.groups()[0])
        if rv == -1.0:
            raise Exception("No output found")
        return rv


def main():
    ML4TADASHI.run(CloverLeaf, {"translator": "Polly"})


def merge01tile(scop_idx, sizex, sizey):
    # merge first two loops and tile
    return [
        [scop_idx, 1, TrEnum.FUSE, 0, 1],
        [scop_idx, 4, TrEnum.FUSE, 0, 1],
        [scop_idx, 3, TrEnum.TILE_2D, sizex, sizey],
    ]


def manual():
    times = [
        0.454686,
        0.282210,
        0.345998,
        0.276447,
        2.708418,
        1.111301,
        2.672818,
        1.125872,
    ]
    indexs = dict([(v, k) for k, v in enumerate(times)])
    max_idx = times.index(max(times))
    # translator = Polly()
    translator = Pet()
    app = CloverLeaf(translator=translator)
    print(f"{len(app.scops)=}")
    app.compile()  # force compile
    otime = 0.0
    print("Measuring: ", end="")
    otime = app.measure(3)
    print(f"{otime=}")
    diffs = [0]
    for base in range(8, 128):
        for dx in diffs:
            for dy in diffs:
                sizex = base + dx
                sizey = sizex + dy
                app.reset_scops()
                trs = []
                for scop_idx in [4, 6]:
                    trs += merge01tile(scop_idx, sizex, sizey)

                    # tile first two loops
                    # trs = [
                    #     [scop_idx, 7, TrEnum.TILE_2D, sizex, sizey],
                    #     [scop_idx, 3, TrEnum.TILE_2D, sizex, sizey],
                    # ]
                app.transform_list(trs)
                # print(app.scops[scop_idx].schedule_tree[0].yaml_str)
                if not app.legal:
                    print("NOT LEGAL")
                    continue
                tapp = app.generate_code()
                ttime = tapp.measure()
                speedup = otime / ttime
                print(f"{(sizex, sizey)=}, {otime=} {ttime=} {speedup=}")


if __name__ == "__main__":
    # main()
    manual()
