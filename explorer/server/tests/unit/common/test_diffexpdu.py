import struct
import unittest

import numpy as np

from server.common.diffexpdu import DiffExArguments, deflate_postings_lists, inflate_postings_lists


class TestDiffexPdu(unittest.TestCase):
    """Test differential expression PDU encoding/decoding."""

    def test_roundtrip_diffex(self):
        """simple round trip"""
        de_args = DiffExArguments(
            mode=DiffExArguments.DiffExMode.TopN,
            params=DiffExArguments.TopNParams(N=50),
            set1=np.arange(0, 100, 2, dtype=np.uint32),
            set2=np.arange(1, 100, 2, dtype=np.uint32),
        )

        encoded = de_args.pack()
        self.assertIsInstance(encoded, bytes)
        self.assertGreater(len(encoded), 0)

        decoded = DiffExArguments.unpack_from(encoded)
        self.assertEqual(decoded, de_args)

    def test_roundtrip_multi_postings_list(self):
        n_obs = 100_000
        n_elem = (50_000, 10, 1000)
        draw = (
            np.random.default_rng().choice(n_obs, sum(n_elem), replace=False).astype(np.uint32)
        )  # uniform distribution

        start = 0
        plists = []
        for n in n_elem:
            arr = draw[start : start + n]
            arr.sort()
            plists.append(arr)
            start = start + n

        buf = deflate_postings_lists(plists, sorted=True)
        rt_plists = inflate_postings_lists(buf)

        for orig, rtrip in zip(plists, rt_plists):
            self.assertTrue(np.array_equal(orig, rtrip))

    def test_roundtrip(self):
        """Test wide variety of densities to trigger all encoding block types"""
        n_obs = 4 * 2**16
        for density in [0.0001, 0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 0.99, 0.9999]:
            n_elem = int(n_obs * density)
            draw = (
                np.random.default_rng().choice(n_obs, n_elem, replace=False).astype(np.uint32)
            )  # uniform distribution
            draw.sort()
            rt = inflate_postings_lists(deflate_postings_lists((draw,), sorted=True))
            self.assertEqual(len(rt), 1)
            self.assertTrue(np.array_equal(draw, rt[0]))

    def test_topN_mode_only(self):
        """Verify that mode check error."""
        with self.assertRaises(AssertionError):
            DiffExArguments(mode=99, params=DiffExArguments.TopNParams(N=1), set1=[0], set2=[1]).pack()

    def test_empty_buffer_error(self):
        """Verify that we catch emtpy buffers."""
        with self.assertRaises(struct.error):
            DiffExArguments.unpack_from(b"")

    def test_error_on_excess_lists(self):
        """Test that an error is raised if number of lists not in range [1, 8]"""
        with self.assertRaises(AssertionError):
            deflate_postings_lists(tuple(), sorted=True)
        with self.assertRaises(AssertionError):
            deflate_postings_lists(([0], [1], [2], [3], [4], [5], [6], [7], [8]), sorted=True)
