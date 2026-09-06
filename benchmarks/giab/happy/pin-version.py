#!/usr/bin/env python
"""Pin archive builds to the checksum-locked hap.py 0.3.15 source version.

The upstream CMake git-describe command produces an empty version outside a git
checkout. Both the C++ and Python generated version files use this one variable.
Compatible with the evaluator's Python 2.7 and host-side Python 3 tests.
"""
from __future__ import print_function

import hashlib
import os
import sys

VERSION = "0.3.15"
EXPECTED_CMAKE_SHA256 = "ba8ef08b590eeb0b0e15756edf17167213db59c7d317795321964ca9b0f115d5"
UPSTREAM_BLOCK = b'''execute_process(COMMAND git describe --tags --always
    OUTPUT_VARIABLE HAPLOTYPES_VERSION
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_STRIP_TRAILING_WHITESPACE
)'''
PINNED_BLOCK = b'''# Rosalind evaluator: version from the checksum-locked upstream release archive.
set(HAPLOTYPES_VERSION "0.3.15")'''


def pinned_cmake(source):
    if hashlib.sha256(source).hexdigest() != EXPECTED_CMAKE_SHA256:
        raise ValueError("hap.py CMakeLists.txt differs from the locked 0.3.15 source")
    if source.count(UPSTREAM_BLOCK) != 1:
        raise ValueError("expected exactly one upstream git-describe version block")
    return source.replace(UPSTREAM_BLOCK, PINNED_BLOCK, 1)


def main():
    if len(sys.argv) != 2:
        raise ValueError("usage: pin-version.py EXTRACTED_HAPPY_SOURCE")
    path = os.path.join(sys.argv[1], "CMakeLists.txt")
    with open(path, "rb") as source:
        patched = pinned_cmake(source.read())
    with open(path, "wb") as output:
        output.write(patched)
    print("Pinned hap.py C++ and Python build version to " + VERSION)


if __name__ == "__main__":
    try:
        main()
    except (IOError, ValueError) as error:
        sys.exit(str(error))
