import math


class Vector(object):
    """
    This class gives simple functionalities related to a Vector
    params:
        coordinates: tuple of vector values
        dimension: Stores the number of dimensions of the vector,
                   equal to the length of coordinates

    Example initialization:
        v = Vector([1, 2, 3])
        v.coordinates => (1, 2, 3)
        v.dimension => 3
    """

    EMPTY_COORDINATES_MSG = "The Coordinates must be nonempty"
    NON_ITERABLE_COORDINATES_MSG = "The Coordinates must be an iterable"
    ZERO_VECTOR_ERROR_MSG = "Cannot normalize zero vector"
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = "No unique parallel component to zero vector"
    NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG = "No unique othogonal component to zero vector"
    ANGLE_WITH_ZERO_VECTOR_MG = "Cannot compute an angle with a zero vector."
    DIMENSION_MORE_THAN_THREE_MSG = "Cross Product can be done on vectors of dimension 3 or less."

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError(Vector.EMPTY_COORDINATES_MSG)
            self.coordinates = tuple(x for x in coordinates)
            self.dimension = len(coordinates)
        except TypeError:
            raise TypeError(Vector.NON_ITERABLE_COORDINATES_MSG)
    
    def __str__(self):
        """
        returns vector as a string:
        Vector (v1, v2, v3,...)
        """
        return "Vector {}".format(self.coordinates)
    
    def __eq__(self, v):
        """
        compares 2 vectors for equality;
        """
        return self.coordinates == v.coordinates


    def __add__(self, v):
        """
        Add two vectors element wise. Addition using + operator
        """
        if self.dimension != v.dimension:
            raise TypeError('Cannot add Vector of dimensions {} and {}'.format(self.dimension, v.dimension))
        return Vector([x + y for x, y in zip(self.coordinates, v.coordinates)])
    
    def __sub__(self, v):
        """
        Subtract two vectors element wise. Subtraction using - operator
        """
        if self.dimension != v.dimension:
            raise TypeError('Cannot subtract Vector of dimensions {} and {}'.format(self.dimension, v.dimension))
        return Vector([x - y for x, y in zip(self.coordinates, v.coordinates)])

    def __mul__(self, n):
        """
        Multiply two vectors element wise. Multiplication using vector1 * scalar
        or vector1 * vector2 operator.
        """
        if isinstance(n, Vector):
            return self.dot(n)
        else:
            return Vector([x * n for x in self.coordinates])

    def __rmul__(self, n):
        """
        Multiply two vectors element wise. Multiplication using scalar * vector1
        or vector2 * vector1 operator.
        """
        if isinstance(n, Vector):
            return self.dot(n)
        else:
            return self.__mul__(n)

    def magnitude(self):
        """
        returns magnitude of a vector
        magnitude = sqrt(v1 * v1 + v2 * v2 + v3 * v3 +...)
        """
        return math.sqrt(sum([x * x for x in self.coordinates]))

    def normalize(self):
        """
        Find unit vector in direction of the current vector
        """
        magnitude = self.magnitude()
        if magnitude == 0:
            raise ZeroDivisionError(Vector.ZERO_VECTOR_ERROR_MSG)
        return Vector([x / magnitude for x in self.coordinates])

    def dot(self, v):
        """
        Find dot product of the two vectors.
        """
        return sum([x * y for x, y in zip(self.coordinates, v.coordinates)])

    def angle_with(self, v, in_degrees=False):
        """
        Find angle between 2 vectors in randians as default.
        if in_degrees = True, then angle returned is in degrees.
        """
        try:
            normal_product = round((self.normalize()).dot(v.normalize()), 5)
            angle_in_rad = math.acos(normal_product)

            if in_degrees:
                return angle_in_rad * 180 / math.pi
            return angle_in_rad
        except Exception as e:
            if str(e) == Vector.ZERO_VECTOR_ERROR_MSG:
                raise Exception(Vector.ANGLE_WITH_ZERO_VECTOR_MG)
            else:
                raise e

    def project(self, b):
        """
        Project current vector on the basis vector b.
        Projection  = magnitude of current vector * unit vector in direction of b.
        """
        try:
            return self.dot(b.normalize()) * b.normalize()
        except Exception as e:
            if str(e) == Vector.ZERO_VECTOR_ERROR_MSG:
                raise Exception(Vector.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def project_orth(self, b):
        """
        Project current vector on a vector perpendicular to the basis vector b.
        current vector = projection on b + projection on orthogonal of b
        => projection on orthogonal of b = current vector - projection on b.
        """
        try:
            parallel_project = self.project(b)
            return self - parallel_project
        except Exception as e:
            if str(e) == Vector.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(Vector.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e

    def __make_3d__(self):
        """
        returns 3D vector (first 3 dimensions) out of current vector
        if dimensions of current vector = 1
        Vector (v1, 0, 0)
        if dimensions of current vector = 2
        Vector (v1, v2, 0)
        if dimensions of current vector >= 3
        Vector (v1, v2, v3)
        """
        if self.dimension > 3:
            return Vector([self.coordinates[i] for i in range(3)])
        li = []
        for i in range(3):
            if i > self.dimension:
                li.append(0)
            else:
                li.append(self.coordinates[i])
        return Vector(li)

    def cross(self, v):
        """
        Find cross product of the first 3D of the current vector with the
        vector v.
        """
        if self.dimension > 3 or v.dimension > 3:
            raise ValueError(Vector.DIMENSION_MORE_THAN_THREE_MSG)
        else:
            print(self)
            print(v)
            v1 = self.__make_3d__()
            v2 = v.__make_3d__()
            x_1, y_1, z_1 = v1.coordinates
            x_2, y_2, z_2 = v2.coordinates

            li = []
            li.append(y_1 * z_2 - y_2 * z_1)
            li.append(-(x_1 * z_2 - z_1 * x_2))
            li.append(x_1 * y_2 - y_1 * x_2)
            return Vector(li)

    def is_parallel_to(self, v):
        """
        Find whether the current vector is parallel to the given vector.
        """
        return (self.is_zero() or
                v.is_zero() or
                self.angle_with(v) == 0 or
                self.angle_with(v) == math.pi)

    def is_orthogonal_to(self, v, tolerance = 1e-10):
        """
        Find whether the current vector is perpendicular to the given vector.
        """
        return abs(self.dot(v)) < tolerance

    def is_zero(self, tolerance=1e-10):
        """
        Find whether the current vector is a zero vector.
        """
        return self.magnitude() < tolerance
    
if __name__ == '__main__':
    my_vector = Vector([1, 2, 3])
    print(my_vector)
    vector1 = Vector([1, 2, 3])
    vector2 = Vector([-1, 2, 3])
    print('Compare vectors: ')
    print(my_vector == vector1)
    print(my_vector == vector2)
    
    print()
    print(my_vector, '+', vector1, '=')
    print(my_vector + vector1)
    
    print()
    print(my_vector, '-', vector1, '=')
    print(my_vector - vector1)
    
    
    print('Magnitude of ', my_vector, '=', my_vector.magnitude())
    print('Unit vector of ', my_vector, '=', my_vector.normalize())
    print('Check for zero vector:', Vector([0, 0, 0]).is_zero())