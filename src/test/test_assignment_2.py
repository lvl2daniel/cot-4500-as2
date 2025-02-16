import unittest
from src.main.assignment_2 import (
    neville_interpolation,
    newton_divided_differences,
    newton_interpolation,
    hermite_interpolation,
    cubic_spline_interpolation
)

class TestInterpolation(unittest.TestCase):
    def test_neville_interpolation(self):
        neville_interpolation([3.6, 3.8, 3.9], [1.675, 1.436, 1.318], 3.7)

    def test_newton_divided_differences(self):
        newton_divided_differences([7.2, 7.4, 7.5, 7.6], [23.5492, 25.3913, 26.8224, 27.4589])

    def test_newton_interpolation(self):
        table = newton_divided_differences([7.2, 7.4, 7.5, 7.6], [23.5492, 25.3913, 26.8224, 27.4589])
        newton_interpolation([7.2, 7.4, 7.5, 7.6], table, 7.3)

    def test_hermite_interpolation(self):
        hermite_interpolation([3.6, 3.8, 3.9], [1.675, 1.436, 1.318], [-1.195, -1.188, -1.182])

    def test_cubic_spline_interpolation(self):
        cubic_spline_interpolation([2, 5, 8, 10], [3, 5, 7, 9])

if __name__ == "__main__":
    unittest.main()
