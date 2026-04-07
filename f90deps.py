#!/usr/bin/env python3
import re
import sys
from collections import defaultdict
from pathlib import Path

# Regexes (case-insensitive)
module_def_re = re.compile(r"^\s*module\s+(\w+)", re.IGNORECASE)
module_proc_re = re.compile(r"^\s*module\s+procedure", re.IGNORECASE)
use_re = re.compile(r"^\s*use\s+(?:,.*::\s*)?(\w+)", re.IGNORECASE)


def find_sources(paths):
    files = []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            files.extend(p.rglob("*.f90"))
        else:
            files.append(p)
    return sorted(set(files))


def scan_file(path):
    provides = set()
    uses = set()

    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()

                # Skip comments
                if line.startswith("!"):
                    continue

                # Ignore "module procedure"
                if module_proc_re.match(line):
                    continue

                m = module_def_re.match(line)
                if m:
                    provides.add(m.group(1).lower())
                    continue

                m = use_re.match(line)
                if m:
                    uses.add(m.group(1).lower())

    except Exception as e:
        print(f"# Warning: failed to read {path}: {e}", file=sys.stderr)

    return provides, uses


def main():
    if len(sys.argv) < 2:
        print("Usage: f90deps.py <files or dirs>")
        sys.exit(1)

    sources = find_sources(sys.argv[1:])

    # Map: module -> file providing it
    module_to_file = {}
    file_provides = {}
    file_uses = {}

    for src in sources:
        provides, uses = scan_file(src)
        file_provides[src] = provides
        file_uses[src] = uses

        for mod in provides:
            module_to_file[mod] = src

    # Generate dependencies
    for src in sources:
        obj = src.with_suffix(".o")
        deps = set()

        for mod in file_uses[src]:
            if mod in module_to_file:
                dep_src = module_to_file[mod]
                if dep_src != src:
                    deps.add(dep_src.with_suffix(".o"))

        if deps:
            dep_list = " ".join(str(d) for d in sorted(deps))
            print(f"{obj}: {dep_list}")


if __name__ == "__main__":
    main()
