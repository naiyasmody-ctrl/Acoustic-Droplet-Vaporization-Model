# Acoustic Droplet Vaporization (ADV) Computational Model

## Overview
Computational physics model predicting when perfluorocarbon (PFP) droplets vaporize under acoustic ultrasound exposure. Built for targeted drug delivery and ultrasound imaging applications.

## Project Summary
This project develops a **quantitative framework** for understanding acoustic droplet vaporization—a critical process in medical ultrasound applications. By modeling how droplet size, temperature, external pressure, and surface coating affect vaporization thresholds, we can optimize acoustic parameters for therapeutic and diagnostic use.

## Key Features
- **54+ simulations** across realistic parameter space
  - Droplet radius: 250 nm - 1000 nm
  - Temperature: 20°C - 37°C (room to body temperature)
  - External pressure: 0.5 - 2.0 atm (tissue depth simulation)
  - Surface tension: 23 - 27 mN/m (coated vs. uncoated)

- **Physics-based predictions** using:
  - Clausius-Clapeyron relation (pressure-temperature coupling)
  - Laplace pressure (surface tension effects)
  - Cavitation mechanics
  - Antoine equation (vapor pressure)

- **Professional visualizations** showing parameter interactions
- **Extensible Python code** ready for Stage 2 (cell binding) and Stage 3 (integrated pipeline)

## Key Findings

### 1. Droplet Size Matters Most
- **Smaller droplets require higher acoustic intensity** due to Laplace pressure (ΔP = 2σ/r)
- 250 nm droplets: ~0.11 MPa threshold
- 500 nm droplets: ~0 MPa (vaporize easily)
- 1000 nm droplets: Spontaneous vaporization at standard conditions

### 2. Body Temperature is Optimal
- **37°C (body temperature) is naturally ideal for PFP vaporization**
- No heating required for in-vivo applications
- Temperature effect validated via Clausius-Clapeyron

### 3. Tissue Depth (Pressure) is the Constraint
- **Higher external pressure increases vaporization threshold**
- Deeper tissue = more acoustic power needed
- Major design consideration for clinical deployment

### 4. Surface Coatings Help
- **Surfactant-coated droplets vaporize more easily** than bare ones
- ~15% reduction in required acoustic intensity
- Improves both stability and efficiency

## Results & Visualizations

### Generated Plots:
1. **threshold_vs_radius.png** - Shows size dependence (Laplace pressure dominance)
2. **threshold_vs_temperature.png** - Shows temperature optimization
3. **threshold_vs_pressure.png** - Shows tissue depth penalty
4. **threshold_3d_surface.png** - Complete parameter landscape

### Data:
- **adv_results.csv** - All 54 simulation results with calculated thresholds

## Repository Structure
