import numpy as np
from scipy.special import erf

def _boys_function_F0(T: np.ndarray) -> np.ndarray:
    """
    Calculate Boys function F0(T)
    F0(T) = 0.5 * sqrt(pi/T) * erf(sqrt(T))

    input:
    T (np.ndarray): scalar or array of T values

    output:
    np.ndarray: F0(T)
    """
    # avoid division by zero
    # F0(T=0) = 1
    sqrt_T = np.sqrt(T)
    T = np.asarray(T, dtype=np.float64)
    if np.any(T < 0):
        raise ValueError("T must be non-negative")
    sqrt_T = np.sqrt(T)

    # threshold below which we use a Taylor expansion for stability
    small_thresh = 1e-8

    # allocate result and compute branches only where needed to avoid warnings/div-by-zero
    result = np.empty_like(T, dtype=np.float64)
    mask_large = T > small_thresh

    if np.any(mask_large):
        result[mask_large] = 0.5 * np.sqrt(np.pi / T[mask_large]) * erf(sqrt_T[mask_large])

    if np.any(~mask_large):
        t = T[~mask_large]
        # Taylor expansion: F0(T) = 1 - T/3 + T^2/10 - T^3/42 + O(T^4)
        result[~mask_large] = 1.0 - t/3.0 + t**2/10.0 - t**3/42.0
    return result

def gaussian_overlap(alpha1: np.float64, alpha2: np.float64, R1: np.ndarray, R2: np.ndarray) -> np.float64:
    """
    Calculate the overlap integral between two Gaussian functions.

    Integration: ∫ N1 * exp(-a1|r-R1|²) * N2 * exp(-a2|r-R2|²) dr

    input:
    alpha1 (float): the exponent parameter of the first Gaussian function
    alpha2 (float): the exponent parameter of the second Gaussian function
    R1 (np.ndarray): the center of the first Gaussian function, shape = (3,)
    R2 (np.ndarray): the center of the first Gaussian function, shape = (3,)

    output:
    float: Overlap
    """
    alpha_p = alpha1 + alpha2
    R1, R2 = np.asarray(R1), np.asarray(R2)
    
    R12_squared = np.sum((R1 - R2)**2)
    
    prefactor = (4 * alpha1 * alpha2 / alpha_p**2)**(3.0 / 4.0)
    exponential = np.exp(-alpha1 * alpha2 / alpha_p * R12_squared)
    
    return prefactor * exponential

## gaussian_coulomb HAS NOT BEEN VERIFIED YET ##
def gaussian_coulomb(alpha1: np.float64, alpha2: np.float64, R1: np.ndarray, R2: np.ndarray, R3: np.ndarray) -> float:
    """
    Calculate the overlap integral between two Gaussian functions and a Coulomb potential. 
    Assumes positive charge.

    Integration: ∫ N1 * exp(-alpha1|r-R1|²) * N2 * exp(-alpha2|r-R2|²) * (1/|r-R3|) / 4pi dr

    input:
    alpha1 (float): the exponent parameter of the first Gaussian function
    alpha2 (float): the exponent parameter of the second Gaussian function
    R1 (np.ndarray): the center of the first Gaussian function, shape = (3,)
    R2 (np.ndarray): the center of the first Gaussian function, shape = (3,)
    R3 (np.ndarray): the position of the charge, shape = (3,)

    output:
    float: Coulomb integral
    """
    alpha_p = alpha1 + alpha2
    R1, R2, R3 = np.asarray(R1), np.asarray(R2), np.asarray(R3)
    
    P = (alpha1 * R1 + alpha2 * R2) / alpha_p
    
    R12_squared = np.sum((R1 - R2)**2)
    PR3_squared = np.sum((P - R3)**2)
    
    T = alpha_p * PR3_squared
    
    prefactor = (4 * alpha1 * alpha2 / np.pi**2)**(3.0 / 4.0) * (2*np.pi/alpha_p)
    exponential = np.exp(-alpha1 * alpha2 / alpha_p * R12_squared)
    
    return prefactor * exponential * _boys_function_F0(T)

# Example usage
if __name__ == "__main__":
    ## example for gaussian_overlap ##
    # STO-3G basis for H and He atoms
    # Format: [alpha, coefficient]
    He_data = [
        [0.6362421394E+01, 0.1543289673E+00],
        [0.1158922999E+01, 0.5353281423E+00],
        [0.3136497915E+00, 0.4446345422E+00]
    ]
    
    H_data = [
        [0.3425250914E+01, 0.1543289673E+00],
        [0.6239137298E+00, 0.5353281423E+00],
        [0.1688554040E+00, 0.4446345422E+00]
    ]
    # Atomic positions (in atomic units)
    # He at origin, H at distance R = 1.4632 a.u. along z-axis
    R_He = np.array([0.0, 0.0, 0.0])
    R_H = np.array([0.0, 0.0, 1.4632])
    
    # Calculate overlap between STO-3G orbitals of He and H
    S_HeH = 0.0
    for alpha_He, coeff_He in He_data:
        for alpha_H, coeff_H in H_data:
            S_HeH += coeff_He * coeff_H * gaussian_overlap(alpha_He, alpha_H, R_He, R_H)
    
    print(f"Overlap between He and H STO-3G orbitals: {S_HeH:.8f}")

    ## another example for gaussian_coulomb ##
    print(f"Coulomb integral between two Gaussians and a charge: {gaussian_coulomb(1, 0.5, np.array([0,0,0]), np.array([0,0,0]), np.array([0,0,0])):.8f}")