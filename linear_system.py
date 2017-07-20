from copy import deepcopy

from vector import Vector
from plane import Plane

class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d
        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def compute_triangular_form(self):
        system = deepcopy(self)

        num_equations = len(system)
        num_variables = self.dimension
        j = 0
        for row in range(num_equations):
            while j < num_variables:
                curr_coef = system[row].normal_vector.coordinates[j]
                if LinearSystem.is_near_zero(curr_coef):
                    is_swapped = system.swap_with_row_below_for_nonzero_coefficient(row, j)
                    if not is_swapped:
                        j += 1
                        continue

                system.clear_coefficients_below(row, j)
                j += 1
                break

        return system

    def swap_with_row_below_for_nonzero_coefficient(self, row_above, coefficient):
        for i in range(row_above, len(self)):
            coef = self[i].normal_vector.coordinates[coefficient]
            if not LinearSystem.is_near_zero(coef):
                self.swap_rows(row_above, i)
                return True
        return False

    def clear_coefficients_below(self, row, coefficient):
        curr_coef = self[row].normal_vector.coordinates[coefficient]
        for i in range(row+1, len(self)):
            row_coef = self[i].normal_vector.coordinates[coefficient]
            if not LinearSystem.is_near_zero(row_coef):
                self.add_multiple_times_row_to_row(-row_coef / curr_coef,
                                                   row, i)

    def compute_rref(self):
        tf = self.compute_triangular_form()

        num_equations = len(self)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()

        for row in range(num_equations-1, -1, -1):
            j = pivot_indices[row]
            if j < 0:
                continue
            tf.scale_row_to_make_coefficient_equal_to_one(row, j)
            tf.clear_coefficients_above(row, j)

        return tf

    def scale_row_to_make_coefficient_equal_to_one(self, row, coefficient):
        coef = self[row].normal_vector.coordinates[coefficient]
        self.multiply_coefficient_and_row(1/coef, row)

    def clear_coefficients_above(self, row, coefficient):
        curr_coef = self[row].normal_vector.coordinates[coefficient]
        for i in range(row-1, -1, -1):
            row_coef = self[i].normal_vector.coordinates[coefficient]
            if not LinearSystem.is_near_zero(row_coef):
                self.add_multiple_times_row_to_row(-row_coef / curr_coef,
                                                   row, i)

    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        n = self.planes[row].normal_vector
        k = self.planes[row].constant_term

        new_normal_vector = n * coefficient
        new_constant_term = k * coefficient

        self.planes[row] = Plane(normal_vector=new_normal_vector, constant_term=new_constant_term)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term

        new_normal_vector = n1 * coefficient + n2
        new_constant_term = k1 * coefficient + k2

        self[row_to_be_added_to] = Plane(normal_vector=new_normal_vector, constant_term=new_constant_term)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector.coordinates)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    @staticmethod
    def is_near_zero(val, eps=1e-10):
        return abs(val) < eps


if __name__ == '__main__':
    p0 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p1 = Plane(normal_vector=Vector([0, 1, 0]), constant_term=2)
    p2 = Plane(normal_vector=Vector([1, 1, -1]), constant_term=3)
    p3 = Plane(normal_vector=Vector([1, 0, -2]), constant_term=2)

    s = LinearSystem([p0, p1, p2, p3])

    print (s)

    print('-'*80)
    print("Testing for 3 helper functions")
    print("swap_rows")
    print("multiply_coefficient_and_row")
    print("add_multiple_times_row_to_row")
    s.swap_rows(0,1)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print('test case 1 failed')

    s.swap_rows(1,3)
    if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
        print('test case 2 failed')

    s.swap_rows(3,1)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print ('test case 3 failed')

    s.multiply_coefficient_and_row(1,0)
    if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
        print ('test case 4 failed')

    s.multiply_coefficient_and_row(-1,2)
    if not (s[0] == p1 and
            s[1] == p0 and
            s[2] == Plane(normal_vector=Vector([-1, -1, 1]), constant_term=-3) and
            s[3] == p3):
        print ('test case 5 failed')

    s.multiply_coefficient_and_row(10,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector([10, 10, 10]), constant_term=10) and
            s[2] == Plane(normal_vector=Vector([-1, -1, 1]), constant_term=-3) and
            s[3] == p3):
        print ('test case 6 failed')

    s.add_multiple_times_row_to_row(0,0,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector([ 10, 10, 10]), constant_term=10) and
            s[2] == Plane(normal_vector=Vector([-1, -1, 1]), constant_term=-3) and
            s[3] == p3):
        print ('test case 7 failed')

    s.add_multiple_times_row_to_row(1,0,1)
    if not (s[0] == p1 and
            s[1] == Plane(normal_vector=Vector([10, 11, 10]), constant_term=12) and
            s[2] == Plane(normal_vector=Vector([-1, -1, 1]), constant_term=-3) and
            s[3] == p3):
        print ('test case 8 failed')

    s.add_multiple_times_row_to_row(-1,1,0)
    if not (s[0] == Plane(normal_vector=Vector([-10, -10, -10]), constant_term=-10) and
            s[1] == Plane(normal_vector=Vector([10, 11, 10]), constant_term=12) and
            s[2] == Plane(normal_vector=Vector([-1, -1, 1]), constant_term=-3) and
            s[3] == p3):
        print ('test case 9 failed')


    print('-'*80)
    print("Testing for compute_triangular_form function")
    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([0, 1, 1]), constant_term=2)
    s = LinearSystem([p1,p2])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2):
        print ('test case 1 failed')

    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=2)
    s = LinearSystem([p1,p2])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == Plane(constant_term=1)):
        print ('test case 2 failed')

    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([0, 1, 0]), constant_term=2)
    p3 = Plane(normal_vector=Vector([1, 1, -1]), constant_term=3)
    p4 = Plane(normal_vector=Vector([1, 0, -2]), constant_term=2)
    s = LinearSystem([p1,p2,p3,p4])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2 and
            t[2] == Plane(normal_vector=Vector([0, 0, -2]), constant_term=2) and
            t[3] == Plane()):
        print ('test case 3 failed')

    p1 = Plane(normal_vector=Vector([0, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([1, -1, 1]), constant_term=2)
    p3 = Plane(normal_vector=Vector([1, 2, -5]), constant_term=3)
    s = LinearSystem([p1,p2,p3])
    t = s.compute_triangular_form()
    if not (t[0] == Plane(normal_vector=Vector([1, -1, 1]), constant_term=2) and
            t[1] == Plane(normal_vector=Vector([0, 1, 1]), constant_term=1) and
            t[2] == Plane(normal_vector=Vector([0, 0, -9]), constant_term=-2)):
        print ('test case 4 failed')


    print('-'*80)
    print("Testing RREF")
    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([0, 1, 1]), constant_term=2)
    s = LinearSystem([p1,p2])
    r = s.compute_rref()
    if not (r[0] == Plane(normal_vector=Vector([1, 0, 0]), constant_term=-1) and
            r[1] == p2):
        print ('test case 1 failed')

    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=2)
    s = LinearSystem([p1,p2])
    r = s.compute_rref()
    if not (r[0] == p1 and
            r[1] == Plane(constant_term=1)):
        print ('test case 2 failed')

    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([0, 1, 0]), constant_term=2)
    p3 = Plane(normal_vector=Vector([1, 1, -1]), constant_term=3)
    p4 = Plane(normal_vector=Vector([1, 0, -2]), constant_term=2)
    s = LinearSystem([p1,p2,p3,p4])
    r = s.compute_rref()
    if not (r[0] == Plane(normal_vector=Vector([1, 0, 0]), constant_term=0) and
            r[1] == p2 and
            r[2] == Plane(normal_vector=Vector([0, 0, -2]), constant_term=2) and
            r[3] == Plane()):
        print ('test case 3 failed')

    p1 = Plane(normal_vector=Vector([0, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([1, -1, 1]), constant_term=2)
    p3 = Plane(normal_vector=Vector([1, 2, -5]), constant_term=3)
    s = LinearSystem([p1,p2,p3])
    r = s.compute_rref()
    if not (r[0] == Plane(normal_vector=Vector([1, 0, 0]), constant_term=23/9) and
            r[1] == Plane(normal_vector=Vector([0, 1, 0]), constant_term=7/9) and
            r[2] == Plane(normal_vector=Vector([0, 0, 1]), constant_term=2/9)):
        print ('test case 4 failed')