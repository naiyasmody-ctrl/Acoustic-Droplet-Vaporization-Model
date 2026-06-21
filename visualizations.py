"""
visualizations.py
Create plots and visualizations from ADV simulation results
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from config import *

def load_results(filename='data/adv_results.csv'):
    """Load results from CSV file."""
    df = pd.read_csv(filename)
    print(f"Loaded {len(df)} results from {filename}")
    return df


def plot_threshold_vs_radius(df, temperature_c=25, static_pressure_atm=1.0):
    """Plot vaporization threshold vs droplet radius."""
    filtered = df[(df['temperature_c'] == temperature_c) & 
                  (df['static_pressure_atm'] == static_pressure_atm)]
    
    grouped = filtered.groupby(['radius_nm', 'surface_tension_mn_m'])['threshold_mpa'].mean().reset_index()
    
    plt.figure(figsize=(10, 6))
    
    for sigma in grouped['surface_tension_mn_m'].unique():
        data = grouped[grouped['surface_tension_mn_m'] == sigma]
        plt.plot(data['radius_nm'], data['threshold_mpa'], 
                marker='o', linewidth=2, markersize=8, 
                label=f'σ = {sigma:.0f} mN/m')
    
    plt.xlabel('Droplet Radius (nm)', fontsize=12, fontweight='bold')
    plt.ylabel('Vaporization Threshold (MPa)', fontsize=12, fontweight='bold')
    plt.title(f'Vaporization Threshold vs Droplet Size\n(T = {temperature_c}°C, P = {static_pressure_atm} atm)', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()
    
    plt.savefig('data/threshold_vs_radius.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: data/threshold_vs_radius.png")
    plt.close()


def plot_threshold_vs_temperature(df, radius_nm=500, static_pressure_atm=1.0):
    """Plot vaporization threshold vs temperature."""
    filtered = df[(df['radius_nm'] == radius_nm) & 
                  (df['static_pressure_atm'] == static_pressure_atm)]
    
    grouped = filtered.groupby(['temperature_c', 'surface_tension_mn_m'])['threshold_mpa'].mean().reset_index()
    
    plt.figure(figsize=(10, 6))
    
    for sigma in grouped['surface_tension_mn_m'].unique():
        data = grouped[grouped['surface_tension_mn_m'] == sigma]
        plt.plot(data['temperature_c'], data['threshold_mpa'], 
                marker='s', linewidth=2, markersize=8, 
                label=f'σ = {sigma:.0f} mN/m')
    
    plt.xlabel('Temperature (°C)', fontsize=12, fontweight='bold')
    plt.ylabel('Vaporization Threshold (MPa)', fontsize=12, fontweight='bold')
    plt.title(f'Vaporization Threshold vs Temperature\n(r = {radius_nm} nm, P = {static_pressure_atm} atm)', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()
    
    plt.savefig('data/threshold_vs_temperature.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: data/threshold_vs_temperature.png")
    plt.close()


def plot_threshold_vs_pressure(df, radius_nm=500, temperature_c=25):
    """Plot vaporization threshold vs static pressure."""
    filtered = df[(df['radius_nm'] == radius_nm) & 
                  (df['temperature_c'] == temperature_c)]
    
    grouped = filtered.groupby(['static_pressure_atm', 'surface_tension_mn_m'])['threshold_mpa'].mean().reset_index()
    
    plt.figure(figsize=(10, 6))
    
    for sigma in grouped['surface_tension_mn_m'].unique():
        data = grouped[grouped['surface_tension_mn_m'] == sigma]
        plt.plot(data['static_pressure_atm'], data['threshold_mpa'], 
                marker='^', linewidth=2, markersize=8, 
                label=f'σ = {sigma:.0f} mN/m')
    
    plt.xlabel('Static Pressure (atm)', fontsize=12, fontweight='bold')
    plt.ylabel('Vaporization Threshold (MPa)', fontsize=12, fontweight='bold')
    plt.title(f'Vaporization Threshold vs Static Pressure\n(r = {radius_nm} nm, T = {temperature_c}°C)', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    plt.tight_layout()
    
    plt.savefig('data/threshold_vs_pressure.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: data/threshold_vs_pressure.png")
    plt.close()


def plot_3d_surface(df):
    """Create a 3D surface plot of threshold across radius and temperature."""
    from mpl_toolkits.mplot3d import Axes3D
    
    filtered = df[(df['static_pressure_atm'] == 1.0) & 
                  (df['surface_tension_mn_m'] == 27)]
    
    pivot = filtered.pivot_table(values='threshold_mpa', 
                                 index='radius_nm', 
                                 columns='temperature_c')
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    X = pivot.columns.values
    Y = pivot.index.values
    X_mesh, Y_mesh = np.meshgrid(X, Y)
    Z = pivot.values
    
    surf = ax.plot_surface(X_mesh, Y_mesh, Z, cmap='viridis', alpha=0.8, edgecolor='none')
    
    ax.set_xlabel('Temperature (°C)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Droplet Radius (nm)', fontsize=11, fontweight='bold')
    ax.set_zlabel('Threshold (MPa)', fontsize=11, fontweight='bold')
    ax.set_title('3D Surface: Vaporization Threshold\n(1 atm, coated droplets)', 
                fontsize=13, fontweight='bold')
    
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    
    plt.tight_layout()
    plt.savefig('data/threshold_3d_surface.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: data/threshold_3d_surface.png")
    plt.close()


def plot_heatmap(df, radius_nm=250):
    """Create a heatmap of threshold across temperature and pressure."""
    import seaborn as sns
    
    filtered = df[(df['radius_nm'] == radius_nm) & 
                  (df['surface_tension_mn_m'] == 27)]
    
    pivot = filtered.pivot_table(values='threshold_mpa', 
                                 index='temperature_c', 
                                 columns='static_pressure_atm')
    
    plt.figure(figsize=(10, 6))
    sns.heatmap(pivot, annot=True, fmt='.4f', cmap='YlOrRd', cbar_kws={'label': 'Threshold (MPa)'})
    
    plt.xlabel('Static Pressure (atm)', fontsize=12, fontweight='bold')
    plt.ylabel('Temperature (°C)', fontsize=12, fontweight='bold')
    plt.title(f'Vaporization Threshold Heatmap\n(r = {radius_nm} nm, coated droplets)', 
              fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    plt.savefig('data/threshold_heatmap.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: data/threshold_heatmap.png")
    plt.close()


def plot_laplace_vs_radius(df):
    """Plot Laplace pressure vs droplet radius."""
    unique = df.drop_duplicates(subset=['radius_nm', 'laplace_pressure_kpa'])
    
    plt.figure(figsize=(10, 6))
    plt.plot(unique['radius_nm'], unique['laplace_pressure_kpa'], 
            marker='o', linewidth=2.5, markersize=10, color='steelblue')
    
    plt.xlabel('Droplet Radius (nm)', fontsize=12, fontweight='bold')
    plt.ylabel('Laplace Pressure (kPa)', fontsize=12, fontweight='bold')
    plt.title('Laplace Pressure vs Droplet Radius\n(ΔP = 2σ/r)', 
              fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plt.savefig('data/laplace_vs_radius.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: data/laplace_vs_radius.png")
    plt.close()


def create_all_visualizations(csv_file='data/adv_results.csv'):
    """Create all visualizations at once."""
    print("\n" + "="*70)
    print("CREATING VISUALIZATIONS")
    print("="*70 + "\n")
    
    df = load_results(csv_file)
    
    print("\n[1/6] Creating threshold vs radius plot...")
    plot_threshold_vs_radius(df, temperature_c=25, static_pressure_atm=1.0)
    
    print("\n[2/6] Creating threshold vs temperature plot...")
    plot_threshold_vs_temperature(df, radius_nm=500, static_pressure_atm=1.0)
    
    print("\n[3/6] Creating threshold vs pressure plot...")
    plot_threshold_vs_pressure(df, radius_nm=500, temperature_c=25)
    
    print("\n[4/6] Creating 3D surface plot...")
    plot_3d_surface(df)
    
    print("\n[5/6] Creating heatmap...")
    plot_heatmap(df, radius_nm=250)
    
    print("\n[6/6] Creating Laplace pressure plot...")
    plot_laplace_vs_radius(df)
    
    print("\n" + "="*70)
    print("ALL VISUALIZATIONS COMPLETE!")
    print("="*70)
    print("\nGenerated files:")
    print("  ✓ data/threshold_vs_radius.png")
    print("  ✓ data/threshold_vs_temperature.png")
    print("  ✓ data/threshold_vs_pressure.png")
    print("  ✓ data/threshold_3d_surface.png")
    print("  ✓ data/threshold_heatmap.png")
    print("  ✓ data/laplace_vs_radius.png")
    print("\n")


if __name__ == "__main__":
    create_all_visualizations()