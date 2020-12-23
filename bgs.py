"""
Recreates line 4 of Table 1 from
Hudson, R. R., and N. L. Kaplan. 1995.
“Deleterious Background Selection with Recombination.”
Genetics 141 (4): 1605–17.

The physical layout of the genome follows their Figure 1.
The fact that the neutral region is non-recombining
is a major simplification for the analysis.  By definition,
a single tree describes the history of that region.  If that
were not the case, we would have to take the weighted
mean of summaries of trees, with tree lengths (relative
to total genome length) as the weights.
"""
import concurrent.futures
import sys

import numpy as np

import freqtracker
import fwdpy11

GENOME_LENGTH = 1.0
R = 0.04
U = 0.08
N = 1600
# number of diploids to sample
NSAM = 10


def runsim(argtuple):
    seed = argtuple

    rng = fwdpy11.GSLrng(seed)

    pdict = {
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0.0, U / 2.0, R),  # The U/2. is from their eqn. 2.
        "nregions": [],
        "sregions": [
            fwdpy11.ConstantS(0, 1.0 / 3.0, 1, -0.02, 1.0),
            fwdpy11.ConstantS(2.0 / 3.0, 1.0, 1, -0.02, 1.0),
        ],
        "recregions": [
            fwdpy11.Region(0, 1.0 / 3.0, 1),
            fwdpy11.Region(2.0 / 3.0, 1.0, 1),
        ],
        "demography": fwdpy11.DiscreteDemography(),
        "simlen": 20 * N,
    }
    params = fwdpy11.ModelParams(**pdict)

    pop = fwdpy11.DiploidPopulation(N, GENOME_LENGTH)

    r = freqtracker.FreqTracker(10 * N)

    fwdpy11.evolvets(
        rng,
        pop,
        params,
        1,
        r,
        track_mutation_counts=True,
        suppress_table_indexing=False,
    )

    return (pop, r)


if __name__ == "__main__":
    seed = int(sys.argv[1])
    pop, r = runsim(seed)

    print(r.trajectories)
