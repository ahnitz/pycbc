""" Tile hierarchy over a parameter space, built online from prior draws.

Each level has a match threshold. A draw descends from the root, compared
at each level only against the current node's children; if the best child
clears that level's threshold the draw descends into it, otherwise the
draw becomes the representative of a new branch, and of a fresh chain of
nodes down the remaining levels. Draws reaching the bottom are recorded against
that leaf's id.

With branching b over d levels a draw costs ~d*b comparisons rather than
N against a flat bank, and each comparison set is one node's children --
small enough to stay cache-resident.

Two metrics, chosen per level
----------------------------
Representatives are stored PSD-whitened and unit-normalized, so

    abs(a.conj() @ b)              is the fixed-time overlap
    max_t abs(IFFT(conj(a) * b))   is the time-maximized match

Coarse levels must use match: a fixed-time overlap badly understates the
similarity of two waveforms that differ mainly by a time shift, so
building with it splits where a real bank would not. At fine thresholds
the two agree closely and overlap is ~19x cheaper.

Because overlap <= match always, `overlap >= threshold` is a sound fast
path at any level -- it implies `match >= threshold` -- so the FFT is
only paid when overlap says "split", which is exactly when it could be
wrong.

Whitening divides by sqrt(PSD) rather than PSD, so that a PAIR of stored
vectors carries exactly one factor of the PSD weighting.
"""
import logging

import numpy


class LevelIndex:
    """ One node's children, with the comparison metric for its level. """

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

    def _prep(self, b, lo, hi):
        A = self.matrix()[lo:hi]
        if A.ndim == 2:                 # (N, nfft) -> single polarization
            A = A[:, None, :]
        return A, numpy.atleast_2d(b), slice(self.kmin, self.kmax)

    def overlap_best(self, b, lo=0, hi=None):
        """ (index, score) of the best fixed-time overlap over reps[lo:hi].
        Score is the MIN over polarizations, so a draw belongs to a tile
        only if every polarization matches its centre. """
        hi = len(self.reps) if hi is None else hi
        if lo >= hi:
            return -1, -1.0
        A, b, sl = self._prep(b, lo, hi)
        ov = numpy.empty((A.shape[0], A.shape[1]))
        for p in range(A.shape[1]):
            ov[:, p] = numpy.abs(A[:, p, sl].conj() @ b[p, sl])
        ov = ov.min(axis=1)
        i = int(numpy.argmax(ov))
        return lo + i, float(ov[i])

    def match_best(self, b, lo=0, hi=None):
        """ Same, for the time-maximized match. One FFT per rep, about 19x
        the cost of the overlap. """
        hi = len(self.reps) if hi is None else hi
        if lo >= hi:
            return -1, -1.0
        A, b, _ = self._prep(b, lo, hi)
        m = numpy.empty((A.shape[0], A.shape[1]))
        for p in range(A.shape[1]):
            Z = numpy.fft.ifft(numpy.conj(A[:, p, :]) * b[p][None, :],
                               axis=1) * self.nfft
            m[:, p] = numpy.abs(Z).max(axis=1)
        m = m.min(axis=1)
        j = int(numpy.argmax(m))
        return lo + j, float(m[j])

    def best(self, b, threshold):
        """ (index, score) of the best child, and whether it clears
        `threshold`. Uses the cheap overlap first; only falls back to the
        FFT-based match when overlap is below threshold AND this level
        wants match, since overlap >= threshold already proves
        match >= threshold.
        """
        i, ov = self.overlap_best(b)
        if ov >= threshold or not self.use_match:
            return i, ov
        return self.match_best(b)


class Node:
    __slots__ = ('rep', 'params', 'children', 'lid')

    def __init__(self, rep, params, lid=-1):
        self.rep = rep
        self.params = params
        self.children = None     # LevelIndex + list of Nodes
        self.lid = lid           # leaf: index into the caller's point blocks


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
        self.n_leaves = 0
        self.n_recheck = 0
        self.n_diverge = 0

    def _new_index(self, depth):
        idx = LevelIndex(self.nfft, self.use_match[depth])
        idx.kmin, idx.kmax = self.kmin, self.kmax
        return idx

    def _grow(self, idx, kids, d, rep, params):
        """ Add a child at level `d`, numbering it if it is a leaf. """
        lid = -1
        if d == len(self.thresholds) - 1:
            lid = self.n_leaves
            self.n_leaves += 1
        # the caller may be handing us a view into a buffer it reuses, and
        # a node keeps its representative for the life of the tree
        rep = numpy.array(rep, copy=True)
        child = Node(rep, params, lid)
        idx.add(rep)
        kids.append(child)
        self.counts[d] += 1
        return child

    def insert(self, rep, params):
        """ Route one draw, creating branches as needed. Returns the id of
        the leaf it landed in.

        Leaves live at a fixed depth, one per threshold, so a node created
        at the last level is a leaf for good and can be numbered on the
        spot. The caller keeps the points themselves in flat arrays keyed
        by that id, which costs far less than a list of tuples per node.
        """
        node = self.root
        for d, thr in enumerate(self.thresholds):
            if node.children is None:
                node.children = (self._new_index(d), [])
            idx, kids = node.children
            if not kids:
                node = self._grow(idx, kids, d, rep, params)
                continue
            i, score = idx.best(rep, thr)
            self.ncomp += len(idx)
            if score < thr:
                node = self._grow(idx, kids, d, rep, params)
            else:
                node = kids[i]
        return node.lid

    def probe(self, rep):
        """ Descend read-only and record what was seen, for a later commit.

        Returns a list with one record per level reached,
        ``(chosen, nkids, i_ov, ov, i_m, m)``. Everything is an index or a
        float, never a node, so a worker can hand it back to the process
        that owns the tree. ``i_m`` is -1 where the match was not needed,
        which is exactly when the overlap already cleared the threshold.

        The list stops at the level that would split, so its length also
        says where that happened.
        """
        node = self.root
        out = []
        for d, thr in enumerate(self.thresholds):
            if node.children is None:
                return out
            idx, kids = node.children
            if not kids:
                return out
            i_ov, ov = idx.overlap_best(rep)
            i_m, m = -1, -1.0
            if ov < thr and idx.use_match:
                i_m, m = idx.match_best(rep)
            if ov >= thr or not idx.use_match:
                chosen, score = i_ov, ov
            else:
                chosen, score = i_m, m
            out.append((chosen, len(kids), i_ov, ov, i_m, m))
            if score < thr:
                return out
            node = kids[chosen]
        return out

    def commit(self, rep, params, probed):
        """ Apply a probed draw, reconciling against anything added since.

        Children are only ever appended, so the ones added after the probe
        are exactly ``kids[nkids:]``. Comparing against just those and
        folding the result into the probed score reproduces what a single
        call over all children would have given: the overlap can only rise,
        so a probe that already cleared the threshold still clears it, and
        where it did not, the match over the union is the better of the two
        halves.

        A record is only reusable while the walk stays on the probed path.
        Diverging invalidates every record below, so the rest is redone in
        full; that costs a normal descent and is rare once the map fills.
        """
        node = self.root
        on_path = True
        for d, thr in enumerate(self.thresholds):
            if node.children is None:
                node.children = (self._new_index(d), [])
            idx, kids = node.children
            if not kids:
                node = self._grow(idx, kids, d, rep, params)
                on_path = False
                continue
            rec = probed[d] if (on_path and d < len(probed)) else None
            if rec is None:
                i, score = idx.best(rep, thr)
                self.ncomp += len(idx)
                if on_path:
                    self.n_diverge += 1
                    on_path = False
            else:
                chosen, nkids, i_ov, ov, i_m, m = rec
                if len(kids) > nkids:
                    self.n_recheck += 1
                    j, ov_new = idx.overlap_best(rep, nkids)
                    if ov_new > ov:
                        i_ov, ov = j, ov_new
                    if ov >= thr or not idx.use_match:
                        i, score = i_ov, ov
                    else:
                        if i_m < 0:      # probe cleared on overlap alone
                            i_m, m = idx.match_best(rep, 0, nkids)
                        j, m_new = idx.match_best(rep, nkids)
                        if m_new > m:
                            i_m, m = j, m_new
                        i, score = i_m, m
                else:
                    i, score = (chosen, ov if (ov >= thr
                                               or not idx.use_match) else m)
            if score < thr:
                node = self._grow(idx, kids, d, rep, params)
                on_path = False
            else:
                if rec is not None and i != rec[0]:
                    on_path = False
                node = kids[i]
        return node.lid

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
