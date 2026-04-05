#!/usr/bin/env python

import os
import re
import sys
from shutil import move, which
from subprocess import CompletedProcess
from typing import Optional

from tadashi.apps import App
from tadashi.translators import Polly, Translator

ml4tadashi = os.path.dirname(__file__)
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
        if src.exists():
            dst = src.parents[1] / src.name
            move(src=src, dst=dst)
        cmd = [
            "make",
            "-j",
            f"SOURCE={self.source.with_suffix('').name}",
            "COMPILER=GNU",
            "C_MPI_COMPILER=mpicc",
            "MPI_COMPILER=mpif90",
        ]
        return cmd

    def run_cmd(self) -> list[str]:
        cmd = [f"./{self.output_binary.name}.x"]
        return cmd

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
    # app = CloverLeaf(translator=Polly("clang"))
    ML4TADASHI.run(CloverLeaf, {"translator": "Polly"})


if __name__ == "__main__":
    main()
