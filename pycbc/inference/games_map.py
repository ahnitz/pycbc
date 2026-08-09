""" Build the tile hierarchy ONLINE from prior draws -- no pregenerated bank.

This replaces the bank.sh -> map.py -> build_tree_reorg.py chain with a
single pass. Each level has a match threshold. A draw descends from the
root, compared at each level ONLY against the current node's children; if
the best child clears that level's threshold it descends, otherwise the
draw becomes the representative of a new branch (and, applying the same
rule down the remaining levels, of a fresh chain to a leaf). Draws
reaching the bottom are stored as that leaf's points.

Cost: with branching b over d levels a draw costs ~d*b comparisons rather
than N against a flat bank, and every comparison set is one node's
children -- small enough to stay cache-resident, which matters more than
the comparison count (a flat 247 MB bank costs 2543 ns/template against
328 ns/template for a cache-resident 400).

Which metric, and why it is not uniform across levels
-----------------------------------------------------
pycbc_brute_bank covers the space using MATCH, maximized over time shift
and phase. map.py instead uses a fixed-time OVERLAP (abs() maximizes
phase only). That was sound where map.py operates -- at minimal_match
0.99 the waveforms are nearly identical and time-maximization buys almost
nothing, so overlap ~= match -- but it fails badly at the coarse levels
this hierarchy introduces. Measured on the narrow config against the real
21-template 0.7 bank, over 300 prior draws:

    time-maximized MATCH : min 0.726, median 0.862, covers 100% at 0.7
    fixed-time OVERLAP   : min 0.047, median 0.669, covers  42% at 0.7

Building with overlap therefore splits where a real bank would not: a
0.7-threshold level came out at 109-123 limbs against the bank's 21, and
ALL 109 had a sibling at match >= 0.7 (median nearest-sibling match
0.891) -- i.e. every limb was redundant under the correct metric.

So: MATCH at coarse levels, OVERLAP at fine levels (where it is a good
approximation and 19x cheaper). Because overlap <= match always,
`overlap >= threshold` is a SOUND fast path at any level -- it implies
match >= threshold -- so the FFT is only paid when overlap says "split",
which is exactly when it could be wrong.
"""
import logging

import numpy


class LevelIndex:
    """ One node's children, with the comparison metric for its level.

    Representatives are stored PSD-whitened and unit-normalized, so
    `abs(a.conj() @ b)` is the fixed-time overlap and
    `max_t |IFFT(conj(a) * b)|` is the time-maximized match. Whitening
    divides by sqrt(PSD) rather than PSD so that a PAIR of stored vectors
    carries exactly one factor of the PSD weighting (dividing by the full
    PSD, as map.py's cached bank rows do, is correct only because the
    other side of that product is a RAW candidate).
    """

    def __init__(self, nfft, use_match):
        self.nfft = int(nfft)
        self.use_match = bool(use_match)
        self.reps = []
        self._mat = None

    def __len__(self):
        return len(self.reps)

    def add(self, rep):
        self.reps.append(rep)
        self._mat = None

    def matrix(self):
        if self._mat is None:
            self._mat = numpy.ascontiguousarray(numpy.array(self.reps))
        return self._mat

    def best(self, b, threshold):
        """ (index, score) of the best child, and whether it clears
        `threshold`. Uses the cheap overlap first; only falls back to the
        FFT-based match when overlap is below threshold AND this level
        wants match, since overlap >= threshold already proves
        match >= threshold.
        """
        A = self.matrix()
        ov = numpy.abs(A[:, self.kmin:self.kmax].conj()
                       @ b[self.kmin:self.kmax])
        i = int(numpy.argmax(ov))
        if ov[i] >= threshold or not self.use_match:
            return i, float(ov[i])
        Z = numpy.fft.ifft(numpy.conj(A) * b[None, :], axis=1) * self.nfft
        m = numpy.abs(Z).max(axis=1)
        j = int(numpy.argmax(m))
        return j, float(m[j])


class Node:
    __slots__ = ('rep', 'params', 'children', 'points')

    def __init__(self, rep, params):
        self.rep = rep
        self.params = params
        self.children = None     # LevelIndex + list of Nodes
        self.points = None       # leaf: list of param tuples


class OnlineTree:
    """ Incrementally-built tile hierarchy.

    Parameters
    ----------
    thresholds : sequence of float
        Match threshold per level, coarsest first (e.g. 0.7, 0.9, 0.99).
    match_below : float
        Levels whose threshold is at or below this use time-maximized
        match; levels above it use the cheaper fixed-time overlap, where
        the two agree closely. Set to 0 to use overlap everywhere, or 1
        to use match everywhere.
    kmin, kmax, nfft : int
        Frequency slice actually occupied by signal, and the padded
        length used for the time-maximization IFFT. nfft = 2*kmax gives
        0.5 ms resolution at delta_f=2, which reproduces
        pycbc.filter.match to ~4e-4.
    """

    def __init__(self, thresholds, kmin, kmax, nfft=None, match_below=0.98):
        self.thresholds = list(thresholds)
        self.kmin, self.kmax = int(kmin), int(kmax)
        self.nfft = int(nfft or 2 * kmax)
        self.use_match = [t <= match_below for t in self.thresholds]
        self.root = Node(None, None)
        self.counts = [0] * len(self.thresholds)
        self.ncomp = 0
        self.n_match_calls = 0
        self.frozen = False

    def _new_index(self, depth):
        idx = LevelIndex(self.nfft, self.use_match[depth])
        idx.kmin, idx.kmax = self.kmin, self.kmax
        return idx

    def insert(self, rep, params):
        """ Route one draw, creating branches as needed. Returns the leaf
        it landed in. """
        node = self.root
        for d, thr in enumerate(self.thresholds):
            if node.children is None:
                node.children = (self._new_index(d), [])
            idx, kids = node.children
            if not kids:
                child = Node(rep, params)
                idx.add(rep)
                kids.append(child)
                self.counts[d] += 1
                node = child
                continue
            i, score = idx.best(rep, thr)
            self.ncomp += len(idx)
            if score < thr:
                child = Node(rep, params)
                idx.add(rep)
                kids.append(child)
                self.counts[d] += 1
                node = child
            else:
                node = kids[i]
        if node.points is None:
            node.points = []
        node.points.append(params)
        return node

    def assign(self, rep):
        """ Route a draw through the FROZEN tree without creating nodes.
        Returns the leaf, or None if the tree cannot accept it (which can
        only happen if the tree is empty). Used for the bulk phase, where
        the structure is fixed and workers are read-only.
        """
        node = self.root
        for d, thr in enumerate(self.thresholds):
            if node.children is None:
                return None
            idx, kids = node.children
            i, _score = idx.best(rep, thr)
            node = kids[i]
        return node

    def leaves(self):
        out = []

        def walk(n):
            if n.children is None:
                out.append(n)
                return
            for c in n.children[1]:
                walk(c)

        for c in (self.root.children[1] if self.root.children else []):
            walk(c)
        return out

    def stats(self):
        b = ['%.1f' % (self.counts[i] / max(self.counts[i - 1], 1))
             for i in range(1, len(self.counts))]
        return {'nodes_per_level': list(self.counts), 'branching': b,
                'leaves': self.counts[-1] if self.counts else 0}


def whiten_padded(h_slice, psd_slice, kmin, kmax, nfft):
    """ PSD-whiten, unit-normalize, zero-pad for the time-maximization
    IFFT. `h_slice` is the raw [kmin:kmax] waveform slice.
    """
    n = abs((h_slice.conj() * h_slice / psd_slice).sum() ** 0.5)
    out = numpy.zeros(nfft, dtype=complex)
    out[kmin:kmax] = (h_slice / numpy.sqrt(psd_slice)) / n
    return out
