#!/usr/bin/env python

import sys
import platform
import argparse
import subprocess
from random import choice, shuffle
from matplotlib.pyplot import Axes
import matplotlib.pyplot as plt
from collections import Counter
from tempfile import TemporaryDirectory
from pathlib import Path
from time import asctime, time

from dark.reads import AARead, DNARead
from dark.aaVars import CODONS
from dark.dimension import dimensionalIterator


def getArgs() -> argparse.Namespace:
    """
    Make an argparse parser and use it to return the command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Compute typical DIAMOND bitscores for a range of identities."
    )

    parser.add_argument(
        "--length",
        default=100,
        type=int,
        metavar="N",
        help="The number of AAs in the test sequences.",
    )

    parser.add_argument(
        "--blastxArgs",
        default="",
        metavar="ARGS",
        help="Additional (non-sensitivity) arguments to pass to 'diamond blastx'.",
    )

    parser.add_argument(
        "--iterations",
        default=10,
        metavar="N",
        type=int,
        help="The number of random sequences to test for each AA identity count.",
    )

    parser.add_argument(
        "--dotsize",
        default=3,
        metavar="N",
        type=int,
        help="The size of the dots for the scatter plots",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Write intermediate processing output to standard error.",
    )

    parser.add_argument(
        "--output",
        metavar="FILENAME",
        default="plot.pdf",
        help="The file to write a plot image to. File format is determined by suffix.",
    )

    parser.add_argument(
        "--errorIncrement",
        default=1,
        metavar="N",
        type=int,
        help="The number of additional errors (non-identical AAs) to add at each step.",
    )

    return parser.parse_args()


def sampleBitScore(
    errorCount: int, length: int, blastxArgs: str, tmpdir: Path, verbose: bool
) -> float:
    """
    Compute a bitscore for two sequences of a given AA length with a given
    number of AA errorCount (mismatches).

    @param errorCount: The C{int} number of AA mismatches.
    @param length: The C{int} length of the AA subject sequence desired.
    @param blastxArgs: Additional arguments to pass to diamond blastx.
    @param tmpdir: The temporary directory in which to operate.
    @param verbose: If C{True} write intermediate processing info to sys.stderr.
    @return: A C{float} bitscore, with a negative value if there is no match.
    """
    aas = list(CODONS)
    subject = AARead("subject", "".join(choice(aas) for _ in range(length)))

    # Make a random list of indices to change.
    indices = list(range(length))
    shuffle(indices)
    indices = indices[:errorCount]

    queryAa = list(subject.sequence)
    for index in indices:
        existingAA = newAA = subject.sequence[index]
        while newAA == existingAA:
            newAA = choice(aas)
        queryAa[index] = newAA

    queryDNA = DNARead("query", "".join(choice(CODONS[aa]) for aa in queryAa))

    # Sanity check that we have the right number of aa mismatches.
    assert errorCount == sum(a != b for (a, b) in zip(queryAa, subject.sequence))

    subjectFile = str(tmpdir / "subject.fasta")
    queryFile = str(tmpdir / "query.fasta")
    dbFile = str(tmpdir / "db")

    # Save the target (subject) sequence.
    with open(subjectFile, "w") as fp:
        print(subject.toString("fasta"), file=fp, end="")

    # And use it to make a DIAMOND database.
    subprocess.check_output(
        f"diamond makedb --in {subjectFile!r} --db {dbFile!r} --quiet", shell=True
    )

    # Save the query sequence.
    with open(queryFile, "w") as fp:
        print(queryDNA.toString("fasta"), file=fp, end="")

    # And match the query against the DIAMOND database.
    match = (
        subprocess.check_output(
            f"diamond blastx {blastxArgs} --query {queryFile!r} "
            f"--db {dbFile!r} --outfmt 6 bitscore",
            shell=True,
            timeout=2,
        )
        .decode("ascii")
        .rstrip()
    )

    if match:
        # Take the first bitscore and check that it's the best.
        bitscores = list(map(float, match.split("\n")))
        bitscore = bitscores.pop(0)
        assert all(bitscore >= b for b in bitscores)
    else:
        # Return zero to indicate no match.
        bitscore = 0.0

    if verbose:
        print("SBJCT:", subject.sequence)
        print("QUERY:", queryAa)
        print("SCORE:", bitscore)
        print("ERROR:", errorCount)

    return bitscore


def plot(
    ax: Axes,
    bottom: bool,
    lhs: bool,
    rhs: bool,
    length: int,
    errorIncrement: int,
    dotsize: int,
    verbose: bool,
    iterations: int,
    blastxArgs: str,
    sensitivity: str,
    tmpdir: Path,
) -> None:
    """
    Make a plot of AA identity vs bitscore for a given sensitivity level.

    @param ax: The C{Axes} in which to plot.
    @param bottom: If C{True}, this is in the bottom row of the subplots,
    @param lhs: If C{True}, this is on the left-hand side of the subplots,
    @param rhs: If C{True}, this is on the right-hand side of the subplots,
    @param length: The length of the amino acid sequences to test.
    @param errorIncrement: The increment to use when increasing the number of
        amino acid mismatches.
    @param dotsize: The size of the scatter plot dots.
    @param verbose: If C{True}, report intermediate progress.
    @param iterations: The number of times to test each non-identical amino
        acid count.
    @param blastxArgs: Extra (non-sensitivity) arguments to pass to DIAMOND.
    @param sensitivity: The DIAMOND sensitivity argument.
    @param tmpdir: The directory in which files for DIAMOND can be written.
    """
    zeroBitscoreScale = 2
    errorCounts = []
    bitscores = []
    color = []
    successCounts: Counter[int] = Counter()
    blastxArgs += f" --{sensitivity}" if sensitivity else ""

    matchColor = "steelblue"
    missColor = "tomato"
    missCounts: Counter[int] = Counter()

    start = time()

    for errorCount in range(0, length, errorIncrement):
        if verbose:
            print(f"Making {errorCount} errors:")

        iteration = 0
        while iteration < iterations:
            if verbose:
                print(f"  Iteration {iteration + 1}/{iterations}")

            try:
                score = sampleBitScore(errorCount, length, blastxArgs, tmpdir, verbose)
            except subprocess.TimeoutExpired:
                # On OS X 14.5 the Python subprocess call to diamond blastx very
                # occasionally does not return. I don't know why. When making the
                # identical call on the command line, diamond (v2.1.9) exits
                # immediately with status zero and no output. So we just do it again.
                print("DIAMOND subprocess timeout! Repeating.", file=sys.stderr)
            else:
                iteration += 1
                if score == 0.0:
                    missCounts[errorCount] += 1
                    bitscores.append(-missCounts[errorCount] * zeroBitscoreScale)
                    color.append(missColor)
                else:
                    bitscores.append(score)
                    color.append(matchColor)
                    successCounts[errorCount] += 1

        errorCounts.extend([errorCount] * iterations)

    elapsed = int(time() - start)

    # Compute the overall success rate for each error count. This will be plotted using
    # the right-hand axis.
    successX = []
    successY = []

    for errorCount in range(0, length, errorIncrement):
        successX.append(errorCount)
        successY.append(successCounts[errorCount] / iterations)

    titleFontSize = 10
    axisFontSize = 8

    title = f"{sensitivity or 'default'} ({elapsed}s)"
    ax.scatter(errorCounts, bitscores, s=dotsize, color=color)
    ax.set_title(title, fontsize=titleFontSize)
    if bottom:
        ax.set_xlabel("AA mismatches", fontsize=axisFontSize)
    if lhs:
        ax.set_ylabel("DIAMOND bitscore", fontsize=axisFontSize)
    ax.set_xlim((0.0, length + 1))
    ax.set_ylim((-(iterations + 1) * zeroBitscoreScale, max(bitscores) + 5.0))
    ax.grid()

    ax2 = ax.twinx()
    ax2.plot(successX, successY, linewidth=1)
    if rhs:
        ax2.set_ylabel("DIAMOND match detection rate", fontsize=axisFontSize)
    else:
        ax2.yaxis.set_visible(False)


def main():
    args = getArgs()

    # The ordering here is according to increasing elapsed time. This is as determined
    # by an earlier run, as opposed to being what a regular English-speaker might
    # expect from the words.
    sensitivities = (
        None,
        "faster",
        "fast",
        "very-sensitive",
        "mid-sensitive",
        "more-sensitive",
        "sensitive",
        "ultra-sensitive",
    )

    cols = 3
    rows = len(sensitivities) // cols + bool(len(sensitivities) % cols)
    fig, axes = plt.subplots(rows, cols, figsize=(12, 9), sharex="col", sharey="row")
    dimensions = dimensionalIterator((rows, cols))

    with TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        for sensitivity in sensitivities:
            print(f"Processing sensitivity: {sensitivity}.")
            row, col = next(dimensions)
            plot(
                axes[row][col],
                row == rows - 1,
                col == 0,
                col == cols - 1,
                args.length,
                args.errorIncrement,
                args.dotsize,
                args.verbose,
                args.iterations,
                args.blastxArgs,
                sensitivity,
                tmpdir,
            )

    # Hide the final subplots (if any) that have no content. We do this because the
    # panel is a rectangular grid and some of the plots at the end of the last row may
    # be unused.
    for row, col in dimensions:
        axes[row][col].axis("off")

    version = (
        subprocess.check_output("diamond --version", shell=True)
        .decode()
        .rstrip()
        .split()[-1]
    )

    fig.subplots_adjust(top=0.87)
    fig.suptitle(
        f"AA mismatch count vs DIAMOND (v{version}) bitscore."
        "\n"
        f"Sequence length {args.length} with {args.iterations} iterations at each "
        "identity level."
        "\n"
        f"Run at {asctime()} on {platform.system()} (release {platform.release()})",
        fontsize=12,
    )
    fig.savefig(args.output)


if __name__ == "__main__":
    main()
