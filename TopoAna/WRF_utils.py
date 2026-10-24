def uv_to_mass(U, V):
    """
    Destagger WRF U and V to mass points.

    Assumes last two dimensions are (y, x).

    U: (..., ny, nx+1)
    V: (..., ny+1, nx)

    Returns
    -------
    Um: (..., ny, nx)
    Vm: (..., ny, nx)
    """

    Um = 0.5 * (U[..., :, :-1] + U[..., :, 1:])
    Vm = 0.5 * (V[..., :-1, :] + V[..., 1:, :])

    return Um, Vm