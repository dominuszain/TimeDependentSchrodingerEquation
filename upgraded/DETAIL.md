# Schrödinger Simulator - Technical Details

## Overview

This code solves the **1D time-dependent Schrödinger equation** for a quantum particle evolving under various potential configurations:

$$i\hbar \frac{\partial \psi}{\partial t} = \left[-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x)\right]\psi$$

The solution uses **atomic units** where ℏ = m = 1, simplifying the equation to:

$$i \frac{\partial \psi}{\partial t} = \left[-\frac{1}{2}\frac{\partial^2}{\partial x^2} + V(x)\right]\psi$$

---

## Numerical Method: Crank-Nicolson

### Why Crank-Nicolson?

The Crank-Nicolson scheme is chosen because it:
1. **Is unconditionally stable** - No CFL-type timestep restriction
2. **Is unitary** - Preserves the norm (probability) exactly
3. **Is second-order accurate** in both space and time
4. **Time-reversible** - Respects quantum mechanics' time-reversal symmetry

### Discretization

The Hamiltonian is discretized using finite differences:

$$\frac{\partial^2 \psi}{\partial x^2} \approx \frac{\psi_{j-1} - 2\psi_j + \psi_{j+1}}{\Delta x^2}$$

The Crank-Nicolson update equation is:

$$\frac{\psi^{n+1} - \psi^n}{\Delta t} = -\frac{i}{\hbar} H \left(\frac{\psi^{n+1} + \psi^n}{2}\right)$$

Rearranging:

$$\left(I + \frac{i\Delta t}{2\hbar} H\right) \psi^{n+1} = \left(I - \frac{i\Delta t}{2\hbar} H\right) \psi^n$$

Define matrices:
- **A** = I + i·dt/(2ℏ) · H  (implicit part)
- **B** = I - i·dt/(2ℏ) · H  (explicit part)

Update rule: **A · ψⁿ⁺¹ = B · ψⁿ**

### Implementation Details

The Hamiltonian matrix H is **tridiagonal** (only nearest-neighbor couplings from the second derivative). This allows:
- Sparse matrix storage (`scipy.sparse.diags`)
- Efficient sparse linear solves (`scipy.sparse.linalg.spsolve`)

The matrices A and B are pre-computed in `_setup_cn_matrices()` for efficiency.

---

## Complex Absorbing Potential (CAP)

### Purpose

Without special treatment, wave packets reflecting off grid boundaries create unphysical interference. The CAP absorbs outgoing waves with minimal reflection.

### Implementation

An imaginary potential is added at both boundaries:

$$V_{CAP}(x) = -i \eta \left(\frac{x}{L}\right)^2$$

where:
- η = CAP strength
- L = CAP width
- The quadratic profile provides smooth turn-on

The CAP is **purely imaginary**, so it removes probability (absorption) without affecting the real energy of absorbed particles.

### Code Location

`_setup_cap()` method in `SchrodingerSolver` class (lines ~140-158)

---

## File Structure

```
schrodinger_sim.py
├── SimulationConfig      # Configuration dataclass
├── SimulationResult      # Results container
├── Potential Functions
│   ├── gaussian_well()
│   └── harmonic_oscillator()
├── Initial State Functions
│   └── gaussian_wave_packet()
├── SchrodingerSolver     # Core numerical solver
│   ├── _setup_cap()
│   ├── _setup_cn_matrices()
│   ├── _normalize()
│   ├── compute_norm()
│   ├── compute_energy_expectation()
│   ├── step()
│   └── run()
└── InteractiveSimulation # Visualization
    ├── PRESETS
    ├── run_preset()
    ├── _create_plot()
    ├── _run_preset_and_close()
    └── _update_plot()
```

---

## Key Variables and Parameters

### SimulationConfig Fields

| Parameter | Default | Description |
|-----------|---------|-------------|
| `hbar` | 1.0 | Reduced Planck constant (atomic units) |
| `mass` | 1.0 | Particle mass (atomic units) |
| `x_min`, `x_max` | 0, 10 | Spatial domain boundaries |
| `num_points` | 300 | Grid resolution |
| `total_time` | 10.0 | Simulation duration |
| `time_step` | 0.001 | Time integration step |
| `cap_enabled` | True | Enable absorbing boundaries |
| `cap_width` | 1.0 | Width of CAP region |
| `cap_strength` | 1.0 | Absorption strength |
| `packet_center` | 4.0 | Initial wave packet position |
| `packet_width` | 0.5 | Initial spatial width (σ) |
| `packet_momentum` | 5.0 | Initial momentum (k₀) |

### Wave Function Representation

The wave function ψ(x) is stored as a complex NumPy array:
- `psi` - Current wave function
- `np.abs(psi)**2` - Probability density |ψ|²
- `np.real(psi)` - Real part
- `np.imag(psi)` - Imaginary part

---

## Physical Interpretation

### What the Simulation Shows

1. **Probability density |ψ|²** - Where the particle is likely to be found
2. **Time evolution** - How the probability distribution changes
3. **Interference patterns** - Quantum superposition effects
4. **Tunneling** - Probability leaking through barriers

### Conservation Laws

- **Norm (total probability)** - Should stay at 1.0; decreases only when CAP absorbs probability
- **Energy expectation value** - Should be approximately constant (small numerical drift possible)

### Diagnostics

The `run()` method prints:
- Final norm (should be ~1.0 if no absorption)
- Energy at start and end (checks conservation)

---

## FAQs

### Q: Why does the probability distribution spread over time?

**A:** This is **quantum wave packet dispersion**. A localized wave packet contains many momentum components (Heisenberg uncertainty). Different momenta travel at different velocities, so the packet spreads. 

The spreading timescale is: t_spread ~ m(Δx)²/ℏ

For electrons: extremely fast. For baseballs: longer than the age of the universe. This is why we don't see it classically.

### Q: Does this spreading mean QM predictions are wrong?

**A:** No - QM is correct! The spreading is real but negligible for macroscopic objects due to:
1. **Large mass** - Spreading time scales with mass
2. **Decoherence** - Environmental interactions destroy quantum coherence

Classical physics emerges as the ℏ → 0 limit (correspondence principle).

### Q: What creates the interference patterns in the Finite Well?

**A:** When the wave packet hits a potential wall, part reflects while the rest continues forward. The **reflected wave interferes with the incoming wave**, creating:
- **Nodes** - Destructive interference (probability ≈ 0)
- **Antinodes** - Constructive interference (enhanced probability)

This is purely quantum - a classical particle would just bounce as a single object.

### Q: Why does the norm sometimes decrease?

**A:** The CAP (Complex Absorbing Potential) removes probability at the boundaries. This is intentional - it simulates an "open" system where particles can escape to infinity. If the packet never reaches the boundaries, norm stays at 1.0.

### Q: How do I know the simulation is physically correct?

**A:** Check these diagnostics:
1. Initial norm ≈ 1.0 (properly normalized)
2. Energy is approximately conserved (constant expectation value)
3. Norm only decreases (never increases) - CAP absorbs, doesn't create
4. Behavior matches physical intuition (oscillation in traps, tunneling through barriers)

### Q: Can I simulate a specific quantum system?

**A:** Yes! Modify the `potential` lambda in a preset or add a new preset. Examples:
- Double well: Two Gaussian wells separated by a barrier
- Step potential: `np.where(x > 5, 10, 0)`
- Periodic potential: `np.sin(2*np.pi*x/period)`

### Q: What if I want to simulate a heavier particle?

**A:** Increase `config.mass`. Heavier particles:
- Spread more slowly
- Tunnel less readily
- Behave more "classically"

### Q: How can I make the packet more particle-like?

**A:** Increase `packet_width`. A wider packet has:
- Smaller momentum uncertainty (Δp ~ ℏ/Δx)
- Slower dispersion
- Sharper "bouncing" behavior

---

## Extending the Code

### Adding a New Preset

```python
"My Custom Potential": {
    "config": SimulationConfig(
        x_min=0.0, x_max=10.0, num_points=300, total_time=8.0,
        time_step=0.001, packet_center=5.0, packet_width=0.5,
        packet_momentum=0.0, cap_enabled=True, cap_width=1.0, cap_strength=0.3,
    ),
    "potential": lambda x: <your potential function>,
    "initial_state": lambda x, cfg: gaussian_wave_packet(x, cfg.packet_center, cfg.packet_width, cfg.packet_momentum),
    "description": "Description of what this shows",
},
```

### Adding a New Potential Function

```python
def my_potential(x: np.ndarray, param: float = 1.0) -> np.ndarray:
    """Description of the potential."""
    return <expression in terms of x>
```

### Exporting Animation

The code uses `FuncAnimation`. To save as MP4:

```python
from matplotlib.animation import FFMpegWriter
writer = FFMpegWriter(fps=30)
anim.save('output.mp4', writer=writer)
```

---

## References

1. Crank, J., & Nicolson, P. (1947). A practical method for solving heat-conduction problems. *Mathematical Proceedings of the Cambridge Philosophical Society*.

2. Griffiths, D. J. (2005). *Introduction to Quantum Mechanics* (2nd ed.). Pearson.

3. Press, W. H., et al. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

---

*Code of Zain Ul Abideen*  
*Upgraded with Claude Code + Qwen3.5*
