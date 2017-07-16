import math
from vector import Vector
from line import Line


class Plane(Line):
    """
    This class gives simple functionalities related to a Plane in 3d
    params:
        normal_vector: of type vector
        constant_term: Constant term in the vector equation the Plane.
        basepoint: generated according to the normal_vector and constant_term
    
    Example of Plane equation:
    Ax + By + Cz = k
    (A, B, C).(x, y) = k
    (A, B, C) => Normal Vector
    k => Constant Term
    
    Example initialization:
    p = Plane(Vector([2, 3, 4]), 5)
    This will generate a plane of the form, 2x + 3y +4z = 5
    """

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'
    
    def __init__(self, normal_vector=None, constant_term=None):
        super(Plane, self).__init__(normal_vector, constant_term)
        self.dimension = 3


if __name__ == '__main__':
    p1 = Plane(Vector([-0.412, 3.806, 0.728]), -3.46)
    p2 = Plane(Vector([1.03, -9.515, -1.82]), 8.65)
    print("Plane \n", p1, "\nand \n", p2)
    print("Parallel: ", p1.is_parallel_to(p2))
    print("Equal: ", p1==p2)
    print('-'*80)

    p1 = Plane(Vector([2.611, 5.528, 0.283]), 4.6)
    p2 = Plane(Vector([7.715, 8.306, 5.342]), 3.76)
    print("Plane \n", p1, "\nand \n", p2)
    print("Parallel: ", p1.is_parallel_to(p2))
    print("Equal: ", p1==p2)
    print('-'*80)

    p1 = Plane(Vector([-7.926, 8.625, -7.212]), -7.952)
    p2 = Plane(Vector([-2.642, 2.875, -2.404]), -2.443)
    print("Plane \n", p1, "\nand \n", p2)
    print("Parallel: ", p1.is_parallel_to(p2))
    print("Equal: ", p1==p2)
    print('-'*80)