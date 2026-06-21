"""
config.py
Configuration and physical constants for ADV model
"""

import numpy as np

# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

# PFP Properties (at 25°C, unless stated)
PFP_BOILING_POINT_ATM = 29.0           # °C at 1 atm
PFP_LIQUID_DENSITY = 1600.0            # kg/m³
PFP_GAS_DENSITY = 0.88                 # kg/m³ at 25°C, 1 atm (approximation)
PFP_HEAT_OF_VAPORIZATION = 22e3        # J/mol (approx 20-25 kJ/mol)
PFP_MOLECULAR_WEIGHT = 0.138           # kg/mol (C5F12)
PFP_SURFACE_TENSION = 0.027            # N/m (27 mN/m, temperature dependent)

# Thermodynamic Constants
GAS_CONSTANT = 8.314                   # J/(mol·K)
ABSOLUTE_ZERO = 273.15                 # K (conversion from Celsius)

# Standard Conditions
STD_PRESSURE = 101325                  # Pa (1 atm)
STD_TEMPERATURE = 25                   # °C
STD_TEMPERATURE_K = STD_TEMPERATURE + ABSOLUTE_ZERO  # Kelvin

# Water/Medium Properties
WATER_DENSITY = 1000.0                 # kg/m³
WATER_VISCOSITY = 0.001                # Pa·s (at 20°C, increases with temp)
WATER_SOUND_SPEED = 1480               # m/s (at 20°C)
WATER_ACOUSTIC_IMPEDANCE = WATER_DENSITY * WATER_SOUND_SPEED  # kg/(m²·s)

# ============================================================================
# DEFAULT SIMULATION PARAMETERS
# ============================================================================

# Droplet Parameters (can be varied in simulations)
DEFAULT_DROPLET_RADIUS = 500e-9        # 500 nm in meters
DEFAULT_SURFACE_TENSION = 0.027        # N/m (coated emulsion)
SURFACE_TENSION_RANGE = [0.020, 0.030] # N/m (uncoated to fully coated)

# Droplet Sizes to Test
DROPLET_RADII_NM = [250, 375, 500, 750, 1000]  # in nm
DROPLET_RADII_M = np.array(DROPLET_RADII_NM) * 1e-9  # converted to meters

# Pressure Parameters
STATIC_PRESSURE_RANGE = [0.5, 1.0, 1.5, 2.0]  # in atm
STATIC_PRESSURE_PA = np.array(STATIC_PRESSURE_RANGE) * STD_PRESSURE  # in Pa

# Temperature Parameters
TEMPERATURE_RANGE = [20, 25, 30, 37, 40]  # °C (room to body temp)

# Acoustic Parameters
ACOUSTIC_FREQUENCY_RANGE = [1e6, 2e6, 3e6, 5e6]  # Hz (1-5 MHz)
ACOUSTIC_PRESSURE_RANGE = np.logspace(4, 6.5, 20)  # 0.01 to ~3 MPa in Pa

# ============================================================================
# ANTOINE EQUATION COEFFICIENTS FOR PFP
# ============================================================================
# Used to calculate vapor pressure as function of temperature
# Formula: log10(P_sat) = A - B / (C + T)
# where P_sat is in bar, T is in °C

ANTOINE_A = 3.9
ANTOINE_B = 870
ANTOINE_C = 235

# ============================================================================
# SIMULATION SETTINGS
# ============================================================================

TOLERANCE = 1e-9                       # Numerical tolerance
MAX_ITERATIONS = 1000                  # For iterative calculations
VERBOSE = True                         # Print debug info
SAVE_RESULTS = True                    # Save to CSV
OUTPUT_DIR = 'data/'                   # Output directory

# Time-domain simulation parameters (if needed later)
TIME_STEP = 1e-6                       # 1 microsecond
MAX_TIME = 100e-6                      # 100 microseconds total simulation
TIME_ARRAY = np.arange(0, MAX_TIME, TIME_STEP)

# ============================================================================
# CONVERSION UTILITIES
# ============================================================================

def pa_to_atm(pressure_pa):
    """Convert pressure from Pa to atm"""
    return pressure_pa / STD_PRESSURE

def atm_to_pa(pressure_atm):
    """Convert pressure from atm to Pa"""
    return pressure_atm * STD_PRESSURE

def celsius_to_kelvin(temp_c):
    """Convert temperature from °C to K"""
    return temp_c + ABSOLUTE_ZERO

def kelvin_to_celsius(temp_k):
    """Convert temperature from K to °C"""
    return temp_k - ABSOLUTE_ZERO

def mpa_to_pa(pressure_mpa):
    """Convert pressure from MPa to Pa"""
    return pressure_mpa * 1e6

def pa_to_mpa(pressure_pa):
    """Convert pressure from Pa to MPa"""
    return pressure_pa / 1e6

# ============================================================================
# VALIDATION
# ============================================================================

if __name__ == "__main__":
    print("ADV Model Configuration Loaded Successfully")
    print(f"PFP Boiling Point at 1 atm: {PFP_BOILING_POINT_ATM}°C")
    print(f"Default Droplet Radius: {DEFAULT_DROPLET_RADIUS*1e9} nm")
    print(f"Default Surface Tension: {DEFAULT_SURFACE_TENSION*1000} mN/m")
    print(f"Water Acoustic Impedance: {WATER_ACOUSTIC_IMPEDANCE} kg/(m²·s)")