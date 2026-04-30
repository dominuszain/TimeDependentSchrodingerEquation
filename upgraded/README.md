# Time-Dependent Schrödinger Equation Solver

A 1D quantum wave packet simulator that solves the time-dependent Schrödinger equation using the Crank-Nicolson finite difference method.

## Features

- **Crank-Nicolson solver** - Unconditionally stable, unitary time evolution
- **Complex Absorbing Potential (CAP)** - Minimizes artificial reflections at boundaries
- **Interactive visualization** - Matplotlib animation with preset scenarios
- **Multiple potential types** - Gaussian wells, harmonic oscillators, finite square wells

## Quick Start

```bash
python schrodinger_sim.py
```

## Presets

1. **Harmonic Trap** - Particle oscillating in a parabolic potential
2. **Gaussian Well** - Particle in a deep Gaussian potential well
3. **Finite Well** - Particle bouncing in a finite square well (shows quantum tunneling and interference)

## Requirements

- Python 3.8+
- NumPy
- Matplotlib
- SciPy

```bash
pip install numpy matplotlib scipy
```

## Usage

The simulation opens an interactive window with buttons to switch between presets. Watch the probability density |ψ|² evolve over time.

---

*Code of Zain Ul Abideen*  
*Upgraded with Claude Code + Qwen3.5*
