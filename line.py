import math
from vector import Vector


class Line(object):
    """
    This class gives simple functionalities related to a line in 2d
    params:
        normal_vector: of type vector
        constant_term: Constant term in the vector equation the line.
        basepoint: generated according to the normal_vector and constant_term
    
    Example of line equation:
    Ax + By = k
    (A, B).(x, y) = k
    (A, B) => Normal Vector
    k => Constant Term
    
    Example initialization:
    l = Line(Vector([2, 3]), 5)
    This will generate a line of the form, 2x + 3y = 5
    """

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'
    
    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2
        if not normal_vector:
            normal_vector = Vector([0]*self.dimension)
        self.normal_vector = normal_vector
        
        if not constant_term:
            constant_term = 0
        self.constant_term = constant_term
        
        self.set_basepoint()

    def set_basepoint(self):
        """
        sets the basepoint of the line using the normal_vector
        and constant_term.

        Use of first non zero coefficient is done to compute the
        basepoint.
        """
        try:
            n = self.normal_vector.coordinates
            c = self.constant_term
            basepoint_coords = [0]*self.dimension

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def find_point_of_intersection(self, l):
        """
        Returns the vector representing the point of
        intersection between two lines.

        returns:
        1. self if the lines are equal
        2. None if the lines are parallel
        3. Else Point of intersection
        """
        try:

            l1 = self.normal_vector.coordinates
            l2 = l.normal_vector.coordinates
            k1 = self.constant_term
            k2 = l.constant_term

            den = l1[0] * l2[1] - l1[1] * l2[0]
            A, B = l1
            C, D = l2
            x_numerator = (D * k1 - B * k2)
            y_numerator = (-C * k1 + A * k2)

            return Vector([x_numerator, y_numerator]) * (1 / den)
        except ZeroDivisionError:
            if self == l:
                return self
            else:
                return None

    def __str__(self):
        """
        Generates the string of equation the line.
        """

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term = False):
            coefficient = round(coefficient, num_decimal_places)

            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'
            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            coefs = n.coordinates
            initial_index = Line.first_nonzero_index(coefs)
            terms = [write_coefficient(coefs[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(coefs[i], num_decimal_places) != 0]
            output = ' '.join(terms)
        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' ={}'.format(constant)

        return output

    def __eq__(self, l):
        """
        returns whether the two lines are equal or not.

        Two lines are equal if the vector from point on one line to
        point on another line is parallel to the normal of any of the line.
        """
        if self.normal_vector.is_zero():
            if not l.normal_vector.is_zero():
                return False
            else:
                diff = self.constant_term - l.constant_term
                # If normal vectors of both lines are zero
                # we check if the constant_term of both lines are eqaul
                # if yes then they are equal.
                return Line.is_near_zero(diff)
        elif l.normal_vector.is_zero():
            return False

        if self.is_parallel_to(l):
            v = self.basepoint - l.basepoint
            return v.is_orthogonal_to(self.normal_vector)
        return False

    def is_parallel_to(self, l):
        """
        returns whether the two lines are parallel or not.

        Two lines are parallel if their normal vectors are parallel.
        """
        return self.normal_vector.is_parallel_to(l.normal_vector)

    @staticmethod
    def is_near_zero(item, tolerance=1e-10):
        """
        Helper method to find that the item is tending to zero or not.
        """
        return abs(item) < tolerance

    @staticmethod
    def first_nonzero_index(iterable):
        """
        Helper method to find first nonzero coefficient index in iterable.
        """
        for k, item in enumerate(iterable):
            if not Line.is_near_zero(item):
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)

if __name__ == '__main__':
    def find_intersection(a, b, k1, c, d, k2):
        l1 = Line(normal_vector = Vector([a, b]), constant_term=k1)
        l2 = Line(normal_vector = Vector([c, d]), constant_term=k2)
        print("intersection of lines: \n", l1, "and\n", l2, "is\n")
        print(l1.find_point_of_intersection(l2))
        print('-'*80)

    find_intersection(4.046, 2.836, 1.21, 10.115, 7.09, 3.025)
    find_intersection(7.204, 3.182, 8.68, 8.172, 4.114, 9.883)
    find_intersection(1.182, 5.562, 6.744, 1.773, 8.343, 9.525)