"""
adv_physics.py
Core physics calculations for acoustic droplet vaporization
"""

import numpy as np
from scipy.optimize import fsolve
from config import *

# ============================================================================
# CLAUSIUS-CLAPEYRON: VAPOR PRESSURE vs TEMPERATURE
# ============================================================================

def vapor_pressure_antoine(temperature_c):
    """
    Calculate saturation vapor pressure of PFP using Antoine equation.
    
    Parameters:
    -----------
    temperature_c : float or array
        Temperature in °C
    
    Returns:
    --------
    float or array
        Saturation vapor pressure in Pa
    """
    # Antoine: log10(P_sat) = A - B / (C + T)
    # Returns P in bar; convert to Pa
    log_p_bar = ANTOINE_A - ANTOINE_B / (ANTOINE_C + temperature_c)
    p_bar = 10 ** log_p_bar
    p_pa = p_bar * 1e5  # 1 bar = 1e5 Pa
    return p_pa


def boiling_point_at_pressure(pressure_pa):
    """
    Find the boiling point (temperature) at a given pressure.
    Uses Antoine equation inverted.
    
    Parameters:
    -----------
    pressure_pa : float
        External pressure in Pa
    
    Returns:
    --------
    float
        Boiling point in °C
    """
    p_bar = pressure_pa / 1e5
    log_p_bar = np.log10(p_bar)
    
    # Rearrange Antoine: T = B / (A - log10(P_sat)) - C
    temperature_c = ANTOINE_B / (ANTOINE_A - log_p_bar) - ANTOINE_C
    return temperature_c


# ============================================================================
# LAPLACE PRESSURE
# ============================================================================

def laplace_pressure(radius_m, surface_tension=DEFAULT_SURFACE_TENSION):
    """
    Calculate Laplace pressure across droplet interface.
    
    Physics:
    --------
    ΔP_Laplace = 2σ/r
    
    This is the pressure difference between inside and outside the droplet.
    It opposes vaporization (makes boiling point higher inside the droplet).
    
    Parameters:
    -----------
    radius_m : float or array
        Droplet radius in meters
    surface_tension : float
        Surface tension in N/m (default: 27 mN/m)
    
    Returns:
    --------
    float or array
        Laplace pressure in Pa
    """
    delta_p = (2 * surface_tension) / radius_m
    return delta_p


# ============================================================================
# VAPORIZATION THRESHOLD
# ============================================================================

def vaporization_threshold_pressure(
    droplet_radius_m,
    static_pressure_pa=STD_PRESSURE,
    temperature_c=STD_TEMPERATURE,
    surface_tension=DEFAULT_SURFACE_TENSION,
    safety_margin=0.0
):
    """
    Calculate the acoustic pressure amplitude required to vaporize a droplet.
    
    Physics Principle:
    ------------------
    For vaporization to occur, the effective local pressure must drop below 
    the vapor pressure of PFP at the given temperature.
    
    During ultrasound rarefaction:
        P_local = P_static - P_acoustic
    
    Vaporization condition:
        P_local < P_vap(T) + ΔP_Laplace
    
    Rearranging:
        P_acoustic > P_static + ΔP_Laplace - P_vap(T) + safety_margin
    
    Parameters:
    -----------
    droplet_radius_m : float or array
        Droplet radius in meters
    static_pressure_pa : float
        Static external pressure in Pa (default: 1 atm)
    temperature_c : float
        Temperature in °C (default: 25°C)
    surface_tension : float
        Surface tension in N/m
    safety_margin : float
        Safety factor to account for nucleation kinetics (default: 0)
        Typical range: 0.05-0.2 MPa
    
    Returns:
    --------
    float or array
        Required acoustic pressure amplitude in Pa
    """
    # Calculate vapor pressure at given temperature
    p_vap = vapor_pressure_antoine(temperature_c)
    
    # Calculate Laplace pressure
    p_laplace = laplace_pressure(droplet_radius_m, surface_tension)
    
    # Vaporization threshold
    p_threshold = static_pressure_pa + p_laplace - p_vap + safety_margin
    
    # Ensure non-negative (shouldn't go negative, but safety check)
    p_threshold = np.maximum(p_threshold, 0)
    
    return p_threshold


def vaporization_threshold_mpa(
    droplet_radius_m,
    static_pressure_pa=STD_PRESSURE,
    temperature_c=STD_TEMPERATURE,
    surface_tension=DEFAULT_SURFACE_TENSION,
    safety_margin=0.0
):
    """
    Same as above but returns result in MPa (more convenient for display)
    """
    p_pa = vaporization_threshold_pressure(
        droplet_radius_m,
        static_pressure_pa,
        temperature_c,
        surface_tension,
        safety_margin
    )
    return p_pa / 1e6


# ============================================================================
# CAVITATION INDEX
# ============================================================================

def cavitation_index(
    static_pressure_pa,
    acoustic_pressure_amplitude_pa,
    vapor_pressure_pa=None,
    temperature_c=STD_TEMPERATURE
):
    """
    Calculate the cavitation index (dimensionless).
    
    Physics:
    --------
    σ_cav = (P_static - P_vapor) / (0.5 * ρ * v²)
    
    Where v is acoustic velocity = P_acoustic / (Z)
    
    Interpretation:
    - Low σ_cav (< 1): High cavitation probability
    - High σ_cav (> 10): Low cavitation probability
    
    Parameters:
    -----------
    static_pressure_pa : float
        Static pressure in Pa
    acoustic_pressure_amplitude_pa : float
        Acoustic pressure amplitude in Pa
    vapor_pressure_pa : float
        Vapor pressure in Pa (if None, calculated from temperature)
    temperature_c : float
        Temperature in °C (used if vapor_pressure_pa is None)
    
    Returns:
    --------
    float
        Cavitation index (dimensionless)
    """
    if vapor_pressure_pa is None:
        vapor_pressure_pa = vapor_pressure_antoine(temperature_c)
    
    # Acoustic velocity
    acoustic_velocity = acoustic_pressure_amplitude_pa / WATER_ACOUSTIC_IMPEDANCE
    
    # Cavitation index
    numerator = static_pressure_pa - vapor_pressure_pa
    denominator = 0.5 * WATER_DENSITY * acoustic_velocity**2
    
    if denominator == 0:
        return np.inf
    
    sigma_cav = numerator / denominator
    return sigma_cav


# ============================================================================
# BUBBLE GROWTH (Post-Vaporization)
# ============================================================================

def bubble_radius_rate_of_change(
    current_radius_m,
    droplet_radius_m,
    acoustic_pressure_pa,
    static_pressure_pa=STD_PRESSURE,
    temperature_c=STD_TEMPERATURE
):
    """
    Calculate the rate of bubble radius change after vaporization.
    
    Simple diffusion-limited growth model:
    dr/dt ∝ (P_vapor - P_inside_bubble)
    
    This is a simplified version. Full treatment requires Stefan problem solution.
    
    Parameters:
    -----------
    current_radius_m : float
        Current bubble radius in meters
    droplet_radius_m : float
        Original droplet radius (for scaling)
    acoustic_pressure_pa : float
        Acoustic pressure amplitude in Pa
    static_pressure_pa : float
        Static pressure in Pa
    temperature_c : float
        Temperature in °C
    
    Returns:
    --------
    float
        dr/dt in m/s
    """
    # Simplified model: driving force is pressure difference
    p_vap = vapor_pressure_antoine(temperature_c)
    p_laplace_bubble = laplace_pressure(current_radius_m, DEFAULT_SURFACE_TENSION)
    
    # Pressure inside bubble (approximate)
    p_inside = static_pressure_pa + p_laplace_bubble
    
    # Driving pressure
    p_drive = p_vap - p_inside
    
    # Growth constant (empirical, depends on diffusivity)
    # Approximate: 10^-3 m/(s·Pa)
    growth_constant = 1e-11
    
    dr_dt = growth_constant * p_drive
    
    return dr_dt


# ============================================================================
# VALIDATION & TESTS
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("ADV Physics Module - Validation Tests")
    print("=" * 70)
    
    # Test 1: Vapor pressure at standard conditions
    print("\n[Test 1] Vapor Pressure at Standard Temperature (25°C)")
    p_vap_25 = vapor_pressure_antoine(25)
    print(f"  P_vap(25°C) = {p_vap_25:.2e} Pa = {p_vap_25/STD_PRESSURE:.3f} atm")
    print(f"  Expected: ~1.15 atm ✓" if 1.1 < p_vap_25/STD_PRESSURE < 1.2 else "  Unexpected value")
    
    # Test 2: Laplace pressure for 500 nm droplet
    print("\n[Test 2] Laplace Pressure for 500 nm Droplet")
    radius_500nm = 500e-9
    p_laplace_500 = laplace_pressure(radius_500nm)
    print(f"  r = 500 nm")
    print(f"  σ = {DEFAULT_SURFACE_TENSION*1000} mN/m")
    print(f"  ΔP_Laplace = {p_laplace_500:.2e} Pa = {p_laplace_500/1e3:.1f} kPa")
    print(f"  Expected: ~100-120 kPa ✓" if 90e3 < p_laplace_500 < 130e3 else "  Unexpected value")
    
    # Test 3: Vaporization threshold
    print("\n[Test 3] Vaporization Threshold for 500 nm Droplet at 1 atm, 25°C")
    p_thresh = vaporization_threshold_mpa(500e-9, STD_PRESSURE, 25)
    print(f"  P_threshold = {p_thresh:.2f} MPa")
    print(f"  Expected: ~0.4-0.6 MPa ✓" if 0.3 < p_thresh < 0.8 else "  Unexpected value")
    
    # Test 4: Threshold variation with droplet size
    print("\n[Test 4] Vaporization Threshold vs Droplet Size (at 1 atm, 25°C)")
    for r_nm in [250, 500, 1000]:
        p_thresh_mpa = vaporization_threshold_mpa(r_nm*1e-9, STD_PRESSURE, 25)
        print(f"  r = {r_nm:4d} nm → P_threshold = {p_thresh_mpa:.3f} MPa")
    
    # Test 5: Cavitation index
    print("\n[Test 5] Cavitation Index at Various Acoustic Pressures")
    for p_acoustic_kpa in [100, 500, 1000]:
        sigma_cav = cavitation_index(STD_PRESSURE, p_acoustic_kpa*1e3)
        print(f"  P_acoustic = {p_acoustic_kpa:4d} kPa → σ_cav = {sigma_cav:.3f}")
    
    print("\n" + "=" * 70)
    print("All tests completed!")
    print("=" * 70)