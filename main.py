"""
main.py
Main script to run Stage 1 ADV simulations
"""

import numpy as np
from adv_model import ADVModel
from adv_physics import (
    vapor_pressure_antoine,
    vaporization_threshold_mpa,
    laplace_pressure
)
from config import *

def main():
    """
    Main execution function for Stage 1 ADV simulations
    """
    
    print("\n" + "=" * 70)
    print("ACOUSTIC DROPLET VAPORIZATION (ADV) MODEL - STAGE 1")
    print("=" * 70)
    
    # Initialize the model
    model = ADVModel(verbose=True)
    
    # ========================================================================
    # PART 1: Compute Full Vaporization Threshold Map
    # ========================================================================
    
    print("\n" + "─" * 70)
    print("PART 1: Full Threshold Map")
    print("─" * 70)
    
    results_df = model.compute_vaporization_threshold_map(
        radii_nm=[250, 500, 1000],
        temperatures_c=[20, 25, 37],
        static_pressures_atm=[0.5, 1.0, 2.0],
        surface_tensions=[0.023, 0.027]  # Uncoated vs. coated
    )
    
    print("\nFirst few results:")
    print(results_df.head(10))
    
    # ========================================================================
    # PART 2: Summary Statistics
    # ========================================================================
    
    print("\n" + "─" * 70)
    print("PART 2: Summary Statistics")
    print("─" * 70)
    
    model.summary_statistics()
    
    # ========================================================================
    # PART 3: Key Insights
    # ========================================================================
    
    print("\n" + "─" * 70)
    print("PART 3: Key Physical Insights")
    print("─" * 70)
    
    # 3.1: Size effect
    print("\n[Insight 1] Droplet Size Effect on Vaporization Threshold")
    print("At 1 atm, 25°C, standard surface tension:\n")
    
    for r_nm in [250, 500, 1000]:
        p_thresh = vaporization_threshold_mpa(r_nm*1e-9, STD_PRESSURE, 25)
        p_laplace = laplace_pressure(r_nm*1e-9, DEFAULT_SURFACE_TENSION) / 1e3
        print(f"  r = {r_nm:4d} nm   →  P_threshold = {p_thresh:.3f} MPa   "
              f"(Laplace: {p_laplace:.1f} kPa)")
    
    print("\nInterpretation: Smaller droplets require HIGHER acoustic pressure")
    print("due to increased Laplace pressure opposing vaporization.")
    
    # 3.2: Temperature effect
    print("\n[Insight 2] Temperature Effect on Vaporization Threshold")
    print("At 1 atm, 500 nm droplet:\n")
    
    for t_c in [20, 25, 37]:
        p_thresh = vaporization_threshold_mpa(500e-9, STD_PRESSURE, t_c)
        p_vap = vapor_pressure_antoine(t_c) / STD_PRESSURE
        print(f"  T = {t_c:2d}°C   →  P_threshold = {p_thresh:.3f} MPa   "
              f"(P_vap = {p_vap:.2f} atm)")
    
    print("\nInterpretation: Higher temperature LOWERS vaporization threshold")
    print("because saturation pressure increases (less acoustic forcing needed).")
    
    # 3.3: Static pressure effect
    print("\n[Insight 3] Static Pressure Effect on Vaporization Threshold")
    print("At 25°C, 500 nm droplet:\n")
    
    for p_atm in [0.5, 1.0, 2.0]:
        p_thresh = vaporization_threshold_mpa(500e-9, atm_to_pa(p_atm), 25)
        print(f"  P_static = {p_atm:.1f} atm   →  P_threshold = {p_thresh:.3f} MPa")
    
    print("\nInterpretation: Higher static pressure INCREASES vaporization threshold")
    print("because it directly opposes vapor formation (Clausius-Clapeyron).")
    
    # ========================================================================
    # PART 4: Save Results
    # ========================================================================
    
    print("\n" + "─" * 70)
    print("PART 4: Saving Results")
    print("─" * 70)
    
    model.save_results()
    print("\n✓ Results saved to data/adv_results.csv")
    
    print("\n" + "=" * 70)
    print("STAGE 1 SIMULATION COMPLETE")
    print("=" * 70)
    print("\nNext Steps:")
    print("  1. Create visualizations.py to plot threshold maps")
    print("  2. Explore acoustic frequency dependence")
    print("  3. Transition to Stage 2: Cell binding dynamics")
    print("\n")


if __name__ == "__main__":
    main()