"""
adv_model.py
Main ADV Model class that organizes Stage 1 simulations
"""

import numpy as np
import pandas as pd
from adv_physics import *
from config import *

class ADVModel:
    """
    Acoustic Droplet Vaporization Model - Stage 1
    
    Simulates vaporization thresholds and dynamics across parameter space.
    """
    
    def __init__(self, verbose=VERBOSE):
        """
        Initialize the ADV model.
        
        Parameters:
        -----------
        verbose : bool
            Print debug information during calculations
        """
        self.verbose = verbose
        self.results = None
        self.threshold_map = None
        
        if self.verbose:
            print("[ADVModel] Initialized")
    
    def compute_vaporization_threshold_map(
        self,
        radii_nm=None,
        temperatures_c=None,
        static_pressures_atm=None,
        surface_tensions=None,
        safety_margin_pa=0
    ):
        """
        Generate a comprehensive map of vaporization thresholds across parameter space.
        
        Parameters:
        -----------
        radii_nm : array-like
            Droplet radii to test in nm (default: config values)
        temperatures_c : array-like
            Temperatures to test in °C (default: config values)
        static_pressures_atm : array-like
            Static pressures to test in atm (default: config values)
        surface_tensions : array-like
            Surface tensions to test in N/m (default: single value)
        safety_margin_pa : float
            Safety margin in Pa to account for nucleation kinetics
        
        Returns:
        --------
        pd.DataFrame
            Table with columns: radius_nm, temperature_c, pressure_atm, 
                              surface_tension, threshold_pa, threshold_mpa
        """
        if radii_nm is None:
            radii_nm = DROPLET_RADII_NM
        if temperatures_c is None:
            temperatures_c = TEMPERATURE_RANGE
        if static_pressures_atm is None:
            static_pressures_atm = STATIC_PRESSURE_RANGE
        if surface_tensions is None:
            surface_tensions = [DEFAULT_SURFACE_TENSION]
        
        results = []
        
        total_iterations = (len(radii_nm) * len(temperatures_c) * 
                          len(static_pressures_atm) * len(surface_tensions))
        iteration = 0
        
        if self.verbose:
            print(f"[ADVModel] Computing threshold map...")
            print(f"  Radii: {len(radii_nm)} values")
            print(f"  Temperatures: {len(temperatures_c)} values")
            print(f"  Static Pressures: {len(static_pressures_atm)} values")
            print(f"  Surface Tensions: {len(surface_tensions)} values")
            print(f"  Total iterations: {total_iterations}")
        
        for r_nm in radii_nm:
            r_m = r_nm * 1e-9
            
            for t_c in temperatures_c:
                
                for p_atm in static_pressures_atm:
                    p_pa = atm_to_pa(p_atm)
                    
                    for sigma in surface_tensions:
                        
                        # Calculate threshold
                        p_threshold_pa = vaporization_threshold_pressure(
                            r_m,
                            p_pa,
                            t_c,
                            sigma,
                            safety_margin_pa
                        )
                        
                        p_threshold_mpa = p_threshold_pa / 1e6
                        
                        results.append({
                            'radius_nm': r_nm,
                            'temperature_c': t_c,
                            'static_pressure_atm': p_atm,
                            'surface_tension_mn_m': sigma * 1000,  # Convert to mN/m
                            'threshold_pa': p_threshold_pa,
                            'threshold_mpa': p_threshold_mpa,
                            'vapor_pressure_atm': vapor_pressure_antoine(t_c) / STD_PRESSURE,
                            'laplace_pressure_kpa': laplace_pressure(r_m, sigma) / 1e3
                        })
                        
                        iteration += 1
                        if self.verbose and iteration % max(1, total_iterations//10) == 0:
                            progress = 100 * iteration / total_iterations
                            print(f"  Progress: {progress:.0f}%")
        
        self.threshold_map = pd.DataFrame(results)
        
        if self.verbose:
            print(f"[ADVModel] Threshold map computed. Shape: {self.threshold_map.shape}")
        
        return self.threshold_map
    
    def save_results(self, filename=None):
        """
        Save results to CSV file.
        
        Parameters:
        -----------
        filename : str
            Output filename (default: data/adv_results.csv)
        """
        if self.threshold_map is None:
            print("[ADVModel] No results to save. Run compute_vaporization_threshold_map() first.")
            return
        
        if filename is None:
            filename = OUTPUT_DIR + 'adv_results.csv'
        
        self.threshold_map.to_csv(filename, index=False)
        
        if self.verbose:
            print(f"[ADVModel] Results saved to {filename}")
    
    def get_threshold_for_conditions(
        self,
        radius_nm,
        temperature_c,
        static_pressure_atm,
        surface_tension=DEFAULT_SURFACE_TENSION
    ):
        """
        Get vaporization threshold for specific conditions (without full map).
        
        Useful for quick lookups.
        
        Parameters:
        -----------
        radius_nm : float
            Droplet radius in nm
        temperature_c : float
            Temperature in °C
        static_pressure_atm : float
            Static pressure in atm
        surface_tension : float
            Surface tension in N/m
        
        Returns:
        --------
        dict
            Dictionary with threshold and other useful info
        """
        r_m = radius_nm * 1e-9
        p_pa = atm_to_pa(static_pressure_atm)
        
        p_threshold = vaporization_threshold_pressure(
            r_m, p_pa, temperature_c, surface_tension
        )
        
        result = {
            'radius_nm': radius_nm,
            'temperature_c': temperature_c,
            'static_pressure_atm': static_pressure_atm,
            'surface_tension_mN_m': surface_tension * 1000,
            'threshold_pa': p_threshold,
            'threshold_mpa': p_threshold / 1e6,
            'threshold_mpa_str': f"{p_threshold/1e6:.2f}",
            'vapor_pressure_atm': vapor_pressure_antoine(temperature_c) / STD_PRESSURE,
            'laplace_pressure_kpa': laplace_pressure(r_m, surface_tension) / 1e3
        }
        
        return result
    
    def summary_statistics(self):
        """
        Print summary statistics of the threshold map.
        """
        if self.threshold_map is None:
            print("[ADVModel] No results computed yet.")
            return
        
        print("\n" + "=" * 70)
        print("THRESHOLD MAP SUMMARY STATISTICS")
        print("=" * 70)
        
        print("\nThreshold (MPa) Statistics:")
        print(self.threshold_map['threshold_mpa'].describe())
        
        print("\n\nThreshold by Droplet Radius:")
        grouped_radius = self.threshold_map.groupby('radius_nm')['threshold_mpa'].agg(['mean', 'min', 'max'])
        print(grouped_radius)
        
        print("\n\nThreshold by Temperature:")
        grouped_temp = self.threshold_map.groupby('temperature_c')['threshold_mpa'].agg(['mean', 'min', 'max'])
        print(grouped_temp)
        
        print("\n\nThreshold by Static Pressure:")
        grouped_pressure = self.threshold_map.groupby('static_pressure_atm')['threshold_mpa'].agg(['mean', 'min', 'max'])
        print(grouped_pressure)
        
        print("\n" + "=" * 70)


if __name__ == "__main__":
    # Quick test of the model
    model = ADVModel(verbose=True)
    
    # Compute full threshold map
    results = model.compute_vaporization_threshold_map()
    
    # Show summary
    model.summary_statistics()
    
    # Save results
    model.save_results()
    
    # Get threshold for specific conditions
    print("\n\nExample: Threshold for 500 nm droplet at 25°C, 1 atm, standard surface tension:")
    specific = model.get_threshold_for_conditions(500, 25, 1.0)
    for key, value in specific.items():
        print(f"  {key}: {value}")