r"""The solution for bicubic interpolation of a 2D function using the following
interpolant extension:

f(x, y) = \sum_{i=0}^{3} \sum_{j=0}^{3} a_{i, j}  x^{i}  y^{j}

Where x and y are the relative coordinates inside a [0, 1] x [0, 1] cell.
"""

import numpy as np

def obtain_the_boundary_condition_coefficient(which:str, x:int, y:int) -> list:
    """Obtain the values of the matrix for different boundary conditions.
    Essentially use the interpolant extension, set the x and y inside and take
    into account which order of derivative we are taking, if we are taking it.
    """
    out = np.zeros(16)
    if which == "f":
        for j in range(4):
            for i in range(4):
                ind = i + 4 * j
                out[ind] = x**(i) * y**(j)
    elif which == "fdx":
        for j in range(4):
            for i in range(4):
                ind = i + 4 * j
                if i == 0:
                    continue
                # In this case we ignore the i == 0
                out[ind] = i * x**(i-1) * y**(j)
    elif which == "fdy":
        for j in range(4):
            for i in range(4):
                ind = i + 4 * j
                if j == 0:
                    continue
                # In this case we ignore the i == 0
                out[ind] = j * x**(i) * y**(j-1)
    elif which == "fdxdy":
        for j in range(4):
            for i in range(4):
                ind = i + 4 * j
                if j == 0 or i == 0:
                    continue
                # In this case we ignore the i == 0
                out[ind] = i * j * x**(i-1) * y**(j - 1)

    return out

def obtain_the_coefficient_matrix():
    """
    With the above defined interpolant extension, what we need is figure out
    the way to calculate the a_{i, j} coefficients. The x and y position inside
    the [0, 1] x [0, 1] cell are determined based on the input data of the
    quantity we wish to interpolate and the point inside the domain for which
    we wish to obtain the interpolated data.

    So in practice we have a rectilinear grid of data, or in other words a
    2D matrix of a quantity we wish to interpolate, F and of course minimally
    two 1D arrays that define the dimensions in the X and Y axis, which define
    the domain area of the quantity.

    Now let's say that we wish to obtain the interpolated values inside a point
    p(X, Y).

    If we have a cell [0, 1] x [0, 1] which is clipped to the cell in our input
    rectiliner grid where the point p(x, y) lies:

                         i   o-------------o
                             |             |
                             |             |
                             |             |
                             |      .      |
                             |      p(X,Y) |
                             |             |
                         i+1 o-------------o
                             j             j+1


    1.) Convert the point location to in-cell position [0, 1] x [0, 1]
    2.) Calculate the coefficients a_{i, j} by solving the system of linear
        equations using boundary conditions in the corner of the cells.
    3.) Obtain the interpolated value at point p.


    1.) and 3.) are straightforward and easy, but in case of 2.) we need to
    solve the system of linear equations:

                            A a = b

    Where b is the vector holding the boundary conditions:

    For boundary conditions we set
                b = [f_{i, j},
                     f_{i+1, j},
                     f_{i, j+1},
                     f_{i+1, j+1},
                     f_dx_{i, j},
                     f_dx_{i+1, j},
                     f_dx_{i, j+1},
                     f_dx_{i+1, j+1},
                     f_dy_{i, j},
                     f_dy_{i+1, j},
                     f_dy_{i, j+1},
                     f_dy_{i+1, j+1},
                     f_dxdy_{i, j},
                     f_dxdy_{i+1, j},
                     f_dxdy_{i, j+1},
                     f_dxdy_{i+1, j+1}]

    The interpolant extensions contains 16 terms of a_{i, j}. Combined with 16
    boundary conditions, we have a matrix A, with size {16, 16}. The boundary
    conditions are in the border of the cell, where the (x, y) is either 1 or
    0.
    """

    # Obtain the A matrix

    import numpy as np
    A = np.zeros((16, 16))

    # First set of boundary conditions
    A[0, :] = obtain_the_boundary_condition_coefficient("f", 0, 0)
    A[1, :] = obtain_the_boundary_condition_coefficient("f", 1, 0)
    A[2, :] = obtain_the_boundary_condition_coefficient("f", 0, 1)
    A[3, :] = obtain_the_boundary_condition_coefficient("f", 1, 1)

    # Second set of boundary conditions - partial derivative dx
    A[4, :] = obtain_the_boundary_condition_coefficient("fdx", 0, 0)
    A[5, :] = obtain_the_boundary_condition_coefficient("fdx", 1, 0)
    A[6, :] = obtain_the_boundary_condition_coefficient("fdx", 0, 1)
    A[7, :] = obtain_the_boundary_condition_coefficient("fdx", 1, 1)

    # Third set of boundary conditions - partial derivative dy
    A[8, :] = obtain_the_boundary_condition_coefficient("fdy", 0, 0)
    A[9, :] = obtain_the_boundary_condition_coefficient("fdy", 1, 0)
    A[10, :] = obtain_the_boundary_condition_coefficient("fdy", 0, 1)
    A[11, :] = obtain_the_boundary_condition_coefficient("fdy", 1, 1)

    # Fourth set of boundary conditions - partial derivative dxdy
    A[12, :] = obtain_the_boundary_condition_coefficient("fdxdy", 0, 0)
    A[13, :] = obtain_the_boundary_condition_coefficient("fdxdy", 1, 0)
    A[14, :] = obtain_the_boundary_condition_coefficient("fdxdy", 0, 1)
    A[15, :] = obtain_the_boundary_condition_coefficient("fdxdy", 1, 1)

    return A

def obtain_the_inverse_matrix(A):
    """Now that we have our matrix A in the A a = b system of linear equations,
    we need to calculate the "inverse" matrix to solve it.
    """

    return np.linalg.inv(A)

def unroll_the_solution(A_inv):
    """If you observe the elements of the A_inv, most of the elements are 0.
    Unrolling it by hand does not yield a lot of code and it is much easier
    for compiler to compile than having for loops.
    """

    for j in range(16):
        out = ""

        for i in range(16):
            v = A_inv[j, i]
            if v == 0:
                continue

            if v < 0:
                # In order to have consistent spacing, make the value
                # positive and prepend the - sign
                out += " - "
                v = -v
            elif out == "":
                # No prefix
                pass
            else:
                out += " + "

            out += f"{v} * b[{i}]"
        print(f"a[{j%4}][{j//4}] = " +out)

    return out

def obtain_expressions_for_function_values_and_derivatives():
    """This function generates the expression for values of the function and
    its derivatives. It goes up to second order.
    """

    # First the value
    print("val = ", end="")
    for j in range(4):
        for i in range(4):
            chunk = ""
            if i != 0:
                chunk += " * " + " * ".join(["x" for i in range(i)])
            if j != 0:
                chunk += " * " + " * ".join(["y" for j in range(j)])

            print(f"a{i + 4*j}{chunk} + /*a[{i}, {j}]*/  \\")

    print()
    print()
    # Now the first derivative dx value
    print("valdx = ", end="")
    for j in range(4):
        for i in range(4):
            chunk = ""
            if i == 0:
                continue
            elif i == 1:
                pass
            else:
                chunk += f" * {i} * " + " * ".join(["x" for _ in range(i-1)])
            if j != 0:
                chunk += " * " + " * ".join(["y" for _ in range(j)])

            print(f"a{i + 4*j}{chunk} + /*a[{i}, {j}]*/  \\")

    print()
    print()
    # Now the first derivative dx value
    print("valdy = ", end="")
    for j in range(4):
        for i in range(4):
            chunk = ""
            if i != 0:
                chunk += " * " + " * ".join(["x" for _ in range(i)])

            if j == 0:
                continue
            elif j == 1:
                pass
            else:
                chunk += f" * {j} * " + " * ".join(["y" for _ in range(j-1)])

            print(f"a{i + 4*j}{chunk} + /*a[{i}, {j}]*/  \\")

    print()
    print()
    # Now the mixed derivative value
    print("valdxdy = ", end="")
    for j in range(4):
        for i in range(4):
            chunk = ""
            if i == 0:
                continue
            elif i == 1:
                pass
            else:
                chunk += f" * {i} * " + " * ".join(["x" for _ in range(i-1)])

            if j == 0:
                continue
            elif j == 1:
                pass
            else:
                chunk += f" * {j} * " + " * ".join(["y" for _ in range(j-1)])

            print(f"a{i + 4*j}{chunk} + /*a[{i}, {j}]*/  \\")

    print()
    print()
    # Now the second derivative dx value
    print("valdxdx = ", end="")
    for j in range(4):
        for i in range(4):
            chunk = ""
            if i <= 1:
                continue
            elif i == 2:
                chunk += f" * 2"
            else:
                chunk += f" * {i*(i-1)} * " + " * ".join(["x" for _ in range(i-2)])

            if j != 0:
                chunk += f" * " + " * ".join(["y" for _ in range(j)])

            print(f"a{i + 4*j}{chunk} + /*a[{i}, {j}]*/  \\")

    print()
    print()
    # Now the second derivative dy value
    print("valdydy = ", end="")
    for j in range(4):
        for i in range(4):
            chunk = ""
            if i != 0:
                chunk += " * " + " * ".join(["x" for _ in range(i)])

            if j <= 1:
                continue
            elif j == 2:
                chunk += f" * 2"
            else:
                chunk += f" * {j*(j-1)} * " + " * ".join(["y" for _ in range(j-2)])

            print(f"a{i + 4*j}{chunk} + /*a[{i}, {j}]*/  \\")

if __name__ == "__main__":

    # The matrix
    A = obtain_the_coefficient_matrix()
    # print(A)

    # The inverse of the matrix.
    A_inv = obtain_the_inverse_matrix(A)

    print(A_inv)

    # Print the unrolled code.
    unroll_the_solution(A_inv)

    obtain_expressions_for_function_values_and_derivatives()