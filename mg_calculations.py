# SM tttt cross section value from MG
sig_SM = 0.009308 * 1000. # fb

# Predefined values of Wilson coefs. for which the EFT tttt cross section was calculated.
# One has to make sure that resulting matrix for sigma_i and sigma_ij is not degenerate
wilson_coefficients = [
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1],
    [1, 1, 1, 1, 1],
    [-1, -1, 1, 1, 1],
    [-1, -1, 1, 0, 1],
    [0, 1, 0, 0, -1],
    [0, 1, 1, 1, 0],
    [1, 0, -1, 1, 0],
    [-1, 0, 0, 1, -1],
    [-1, 0, 0, -1, 1],
    [0, 1, -1, 1, -1],
    [0, 1, 0, -1, 0],
    [0, 0, -1, -1, 1],
    [1, -1, 0, -1, 0],
    [1, 1, 0, -1, 1],
    [0, 1, 0, -1, 1],
    [1, -1, -1, 0, 1]
]

# EFT cross section values for different values of Ci in the vectors c1, c2, c3, ... c19, c20
MG_SM = [0.01557, 0.01564, 0.0102, 0.01116, 0.01022, 0.02873, 0.0203, 0.01704, 0.01527, 0.02066, 0.0168, 0.01741,
         0.01567, 0.01342, 0.01839, 0.01208, 0.02113, 0.0283, 0.01983, 0.02386]  # pb
