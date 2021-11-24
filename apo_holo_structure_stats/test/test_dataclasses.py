from unittest import TestCase

from apo_holo_structure_stats.core.dataclasses import contained_in_segment


class TestDataclasses(TestCase):
    def test_contained_in_segment(self):
        segment_starts = [20, 50, 60, 120]
        segment_ends =   [30, 55, 119, 150]

        # edges
        self.assertFalse(contained_in_segment(1, segment_starts, segment_ends))
        self.assertFalse(contained_in_segment(151, segment_starts, segment_ends))

        self.assertTrue(contained_in_segment(20, segment_starts, segment_ends))
        self.assertTrue(contained_in_segment(30, segment_starts, segment_ends))
        self.assertTrue(contained_in_segment(60, segment_starts, segment_ends))
        self.assertTrue(contained_in_segment(119, segment_starts, segment_ends))
        self.assertTrue(contained_in_segment(120, segment_starts, segment_ends))
        self.assertTrue(contained_in_segment(150, segment_starts, segment_ends))

        # inner
        self.assertFalse(contained_in_segment(49, segment_starts, segment_ends))
        self.assertFalse(contained_in_segment(31, segment_starts, segment_ends))

        self.assertTrue(contained_in_segment(52, segment_starts, segment_ends))
        self.assertTrue(contained_in_segment(52, segment_starts, segment_ends))

