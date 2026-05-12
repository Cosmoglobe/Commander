from __future__ import annotations
import numpy as np

P = [
    -0.27292696,
    -0.07629969,
    -0.02819452,
    -0.22797056,
    -0.01471565,
    0.27058160,
    0.54852384,
    0.48051509,
    -0.56800938,
    -0.60441560,
    -0.62930065,
    -1.74114454,
    0.30803317,
    1.50880086,
    0.93412077,
    0.25795794,
    1.71547508,
    0.98938102,
    -0.93678576,
    -1.41601920,
    -0.63915306,
    0.02584375,
    -0.53022337,
    -0.83180469,
    0.08693841,
    0.33887446,
    0.52032238,
    0.14381585,
]


IT28 = 2**28
N_10 = 1024


def deproject_face(x: float, y: float) -> float:
    x_squared = x**2
    y_squared = y**2

    return x * (
        1
        + (1 - x_squared)
        * (
            P[0]
            + x_squared
            * (
                P[1]
                + x_squared
                * (
                    P[3]
                    + x_squared
                    * (
                        P[6]
                        + x_squared * (P[10] + x_squared * (P[15] + x_squared * P[21]))
                    )
                )
            )
            + y_squared
            * (
                P[2]
                + x_squared
                * (
                    P[4]
                    + x_squared
                    * (
                        P[7]
                        + x_squared * (P[11] + x_squared * (P[16] + x_squared * P[22]))
                    )
                )
                + y_squared
                * (
                    P[5]
                    + x_squared
                    * (
                        P[8]
                        + x_squared * (P[12] + x_squared * (P[17] + x_squared * P[23]))
                    )
                    + y_squared
                    * (
                        P[9]
                        + x_squared * (P[13] + x_squared * (P[18] + x_squared * P[24]))
                        + y_squared
                        * (
                            P[14]
                            + x_squared * (P[19] + x_squared * P[25])
                            + y_squared
                            * (P[20] + x_squared * P[26] + y_squared * P[27])
                        )
                    )
                )
            )
        )
    )


def get_tangent_plane_coords(x: float, y: float) -> tuple[float, float]:
    """
    Args:
        x: database coordinate in range -1 to +1
        y: database coordinate in range -1 to +1

    Returns:
        x_i: tangent plane coordinate in range -1 to +1
        eta: tangent plane coordinate in range -1 to +1
    """
    x_i = deproject_face(x, y)
    eta = deproject_face(y, x)

    return x_i, eta


def cube_face_to_vector(
    n_face: int, x_i: float, eta: float
) -> tuple[float, float, float]:
    # ! CONVERTS FACE NUMBER n_f (0-5) AND  x_i, eta: float (-1. - +1.)
    # ! INTO A UNIT VECTOR C

    # ! Adapted from xyaxis.f which was extracted from COBE's upx_pixel_vector.f
    # ! Converted to Fortran 90, August 1998.
    # ! A.J. Banday, MPA.

    # ! To preserve symmetry, the normalization sum must
    # ! always have the same ordering (i.e. largest to smallest).

    x_i_1 = max(x_i_abs := abs(x_i), eta_abs := abs(eta))
    eta_1 = min(x_i_abs, eta_abs)
    norm = 1 / np.sqrt(1 + x_i_1**2 + eta_1**2)

    idx = n_face + 1
    if idx == 1:
        return (norm, -eta * norm, x_i * norm)
    elif idx == 2:
        return (norm, eta * norm, x_i * norm)
    elif idx == 3:
        return (-x_i * norm, norm, eta * norm)
    elif idx == 4:
        return (-x_i * norm, -norm, eta * norm)
    elif idx == 5:
        return (x_i * norm, -norm, eta * norm)
    elif idx == 6:
        return (x_i * norm, norm, -eta * norm)

    raise ValueError(f"Invalid face number: {n_face}")


def pix2vec(ipix: int) -> tuple[float, float, float]:
    """Converts a res 15 quadcube pixel to a unit vector."""

    i_x = np.zeros(N_10)
    i_y = np.zeros(N_10)

    for kpix in range(N_10):
        jpix = kpix
        ix = 0
        iy = 0
        ip = 1
        while jpix != 0:
            id = jpix % 2
            jpix = jpix // 2
            ix = id * ip + ix
            id = jpix % 2
            jpix = jpix // 2
            iy = id * ip + iy
            ip = 2 * ip

        i_x[kpix] = ix
        i_y[kpix] = iy

    n_face = ipix // IT28
    n = ipix % IT28
    i = n % N_10
    n = n // N_10
    j = n % N_10
    k = n // N_10
    j_x = N_10 * i_x[k] + 32 * i_x[j] + i_x[i]
    y_y = N_10 * i_y[k] + 32 * i_y[j] + i_y[i]
    x = (j_x - 8191.5) / 8192.0  # Database coordinates.  Range
    y = (y_y - 8191.5) / 8192.0  # -1 -> 1 covers the square face

    x_i, eta = get_tangent_plane_coords(x, y)  # Distort to tangent plane
    # TWO ITERATIONS TO IMPROVE THE DISTORTION
    for _ in range(2):
        x_p, y_p = get_tangent_plane_coords(x_i, eta)
        x_i -= x_p - x
        eta -= y_p - y

    return cube_face_to_vector(n_face, x_i, eta)
