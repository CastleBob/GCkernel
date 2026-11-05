"""Read-only exact arithmetic audit of SL0's captured transformed plane values.

Run with Python 3; standard library only. This is NOT a kernel acceptance test
or an alternative Boolean implementation. Float values are treated as exact
binary inputs using Fraction, so this cannot recover pre-transform intent.
Captured by temporary logging in endcut::hit at kernel 97f06f4 while running
the existing channel_coincidence Rust test. Logging was removed afterward.
"""
from fractions import Fraction
import unittest


# (px, py, nx, ny), in the transverse frame at the refused band.
PLANES = (
    (1.2532904107566987, 0.364967378145554,
     -0.27959354996155233, 0.9601184545772968),
    (-1.5249931406297772e-15, -4.440892098500626e-16,
     0.2795935499615516, -0.9601184545772972),
)
LO = Fraction(9.53120712893607931e-16)
HI = Fraction(0.331601948361131105)


def exact_coefficients(plane):
    px, py, nx, ny = map(Fraction, plane)
    return -nx / ny, py + nx * px / ny


def stored_coefficients(plane):
    px, py, nx, ny = plane
    return -nx / ny, py + nx * px / ny


class SourcePlaneAudit(unittest.TestCase):
    def test_capture_reproduces_logged_graph_values(self):
        self.assertEqual(stored_coefficients(PLANES[0]),
                         (0.2912073490813658, -2.7755575615628914e-16))
        self.assertEqual(stored_coefficients(PLANES[1]),
                         (0.291207349081365, 0.0))

    def test_exact_stored_supports_cross_inside_refused_band(self):
        am, ab = exact_coefficients(PLANES[0])
        bm, bb = exact_coefficients(PLANES[1])
        self.assertNotEqual(am, bm)
        self.assertNotEqual(ab, bb)
        self.assertLess((am - bm) * LO + ab - bb, 0)
        self.assertGreater((am - bm) * HI + ab - bb, 0)
        intersection = (bb - ab) / (am - bm)
        self.assertLess(LO, intersection)
        self.assertLess(intersection, HI)

    def test_conversion_error_is_small_fraction_of_source_intercept_gap(self):
        exact = tuple(map(exact_coefficients, PLANES))
        stored = tuple(map(stored_coefficients, PLANES))
        source_gap = abs(exact[0][1] - exact[1][1])
        conversion_error = sum(abs(Fraction(s[1]) - e[1])
                               for s, e in zip(stored, exact))
        # Diagnostic ratio, not an acceptance tolerance: most of the gap was
        # already in the input planes before this local conversion.
        self.assertLess(conversion_error / source_gap, Fraction(7, 100))


if __name__ == "__main__":
    unittest.main()
