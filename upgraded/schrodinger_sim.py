"""
Time-dependent Schrödinger equation solver in 1D.

Uses Crank-Nicolson finite difference method with complex absorbing potential (CAP)
boundary conditions for numerical stability and minimal artificial reflections.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from dataclasses import dataclass
from typing import Callable, Optional
import warnings
from matplotlib.widgets import Button

# Handle numpy integration function naming (changed across versions)
_numeric_integral = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', np.trapezoid)


@dataclass
class SimulationConfig:
    """Configuration parameters for the Schrödinger simulation."""
    # Physical constants (atomic units: ħ = m = 1)
    hbar: float = 1.0
    mass: float = 1.0

    # Grid parameters
    x_min: float = 0.0
    x_max: float = 10.0
    num_points: int = 300
    total_time: float = 10.0

    # Time stepping
    time_step: float = 0.001

    # Absorbing boundary (CAP) parameters
    cap_enabled: bool = True
    cap_width: float = 1.0  # Width of absorbing region at each boundary
    cap_strength: float = 1.0  # Strength of absorbing potential

    # Initial wave packet parameters
    packet_center: float = 4.0
    packet_width: float = 0.5
    packet_momentum: float = 5.0  # Initial momentum (positive = moving right)


@dataclass
class SimulationResult:
    """Container for simulation results and diagnostics."""
    x: np.ndarray
    time_axis: np.ndarray
    probability_density: np.ndarray
    real_part: np.ndarray
    imag_part: np.ndarray
    norm_history: np.ndarray
    energy_history: np.ndarray


# =============================================================================
# Potential Functions
# =============================================================================

def gaussian_well(x: np.ndarray, depth: float = 10.0, center: float = 5.0, width: float = 2.0) -> np.ndarray:
    """Create a Gaussian potential well."""
    return -depth * np.exp(-((x - center) ** 2) / (2 * width ** 2))


def harmonic_oscillator(x: np.ndarray, k: float = 1.0, center: float = 5.0) -> np.ndarray:
    """Create a harmonic oscillator potential."""
    return 0.5 * k * (x - center) ** 2


# =============================================================================
# Initial Wave Functions
# =============================================================================

def gaussian_wave_packet(x: np.ndarray, center: float = 4.0,
                         width: float = 0.5, momentum: float = 5.0) -> np.ndarray:
    """
    Create a Gaussian wave packet with initial momentum.

    ψ(x) = exp(-(x-x₀)²/(2σ²)) * exp(i*k₀*x)
    """
    envelope = np.exp(-((x - center) ** 2) / (2 * width ** 2))
    phase = np.exp(1j * momentum * x)
    return envelope * phase


# =============================================================================
# Core Solver
# =============================================================================

class SchrodingerSolver:
    """
    Solves the 1D time-dependent Schrödinger equation using Crank-Nicolson method.

    The Crank-Nicolson scheme is unconditionally stable and unitary (preserves norm),
    making it ideal for quantum time evolution.
    """

    def __init__(self, config: SimulationConfig,
                 potential_func: Callable[[np.ndarray], np.ndarray] = None,
                 initial_state_func: Callable[[np.ndarray], np.ndarray] = None):
        self.config = config
        self.x = np.linspace(config.x_min, config.x_max, config.num_points)
        self.dx = self.x[1] - self.x[0]

        # Set up potential
        if potential_func is None:
            self.V = gaussian_well(self.x).astype(complex)
        else:
            self.V = potential_func(self.x).astype(complex)

        # Set up absorbing boundary potential
        self._setup_cap()

        # Set up initial state
        if initial_state_func is None:
            self.psi = gaussian_wave_packet(
                self.x,
                config.packet_center,
                config.packet_width,
                config.packet_momentum
            )
        else:
            self.psi = initial_state_func(self.x)

        # Normalize initial state
        self.psi = self._normalize(self.psi)

        # Pre-compute Crank-Nicolson matrices
        self._setup_cn_matrices()

        # Diagnostics
        self.norm_history = []
        self.energy_history = []

    def _setup_cap(self):
        """Set up complex absorbing potential at boundaries."""
        self.cap = np.zeros_like(self.x, dtype=complex)

        if not self.config.cap_enabled:
            return

        x = self.x
        L = self.config.cap_width
        eta = self.config.cap_strength

        # Left boundary absorbing region (increases from 0 at x=L to max at x=0)
        left_mask = x < L
        self.cap[left_mask] = -1j * eta * ((L - x[left_mask]) / L) ** 2

        # Right boundary absorbing region (increases from 0 at x=x_max-L to max at x=x_max)
        right_boundary = self.x[-1] - L
        right_mask = x >= right_boundary
        self.cap[right_mask] = -1j * eta * ((x[right_mask] - right_boundary) / L) ** 2

    def _setup_cn_matrices(self):
        """
        Pre-compute the Crank-Nicolson matrices A and B.

        The update equation is: A @ psi^{n+1} = B @ psi^n
        where A = I + i*dt/(2ħ) * H and B = I - i*dt/(2ħ) * H
        """
        dt = self.config.time_step
        hbar = self.config.hbar
        m = self.config.mass
        dx = self.dx

        # Kinetic energy prefactor
        alpha = 1j * dt * hbar / (4 * m * dx ** 2)
        beta = 1j * dt / (2 * hbar)

        n = len(self.x)

        # Build Hamiltonian matrix (sparse tridiagonal)
        # H = -ħ²/(2m) * d²/dx² + V(x) + CAP
        diag_V = self.V + np.real(self.cap)  # Real part of potential
        imag_cap = np.imag(self.cap)  # Imaginary absorbing potential

        # Matrix A (implicit part)
        A_diag = np.ones(n) + alpha * 2 + beta * (diag_V + 1j * imag_cap)
        A_offdiag = -alpha * np.ones(n - 1)

        # Matrix B (explicit part)
        B_diag = np.ones(n) - alpha * 2 - beta * (diag_V + 1j * imag_cap)
        B_offdiag = alpha * np.ones(n - 1)

        # Create sparse tridiagonal matrices
        from scipy.sparse import diags
        from scipy.sparse.linalg import spsolve

        self.A = diags([A_offdiag, A_diag, A_offdiag], [-1, 0, 1], format='csc')
        self.B = diags([B_offdiag, B_diag, B_offdiag], [-1, 0, 1], format='csc')
        self.spsolve = spsolve

    def _normalize(self, psi: np.ndarray) -> np.ndarray:
        """Normalize a wave function."""
        norm = np.sqrt(_numeric_integral(np.abs(psi) ** 2, self.x))
        if norm > 0:
            return psi / norm
        return psi

    def compute_norm(self, psi: np.ndarray) -> float:
        """Compute the norm of a wave function using numerical integration."""
        return _numeric_integral(np.abs(psi) ** 2, self.x)

    def compute_energy_expectation(self, psi: np.ndarray) -> float:
        """Compute the expectation value of energy."""
        hbar = self.config.hbar
        m = self.config.mass
        dx = self.dx

        # Kinetic energy: -ħ²/(2m) * d²ψ/dx²
        d2psi = np.zeros_like(psi)
        d2psi[1:-1] = (psi[:-2] - 2 * psi[1:-1] + psi[2:]) / dx ** 2
        # 2nd order one-sided BC (safe for small grids)
        d2psi[0] = (2 * psi[0] - 5 * psi[1] + 4 * psi[2] - psi[3]) / dx ** 2 if len(psi) > 3 else 0
        d2psi[-1] = (2 * psi[-1] - 5 * psi[-2] + 4 * psi[-3] - psi[-4]) / dx ** 2 if len(psi) > 3 else 0

        kinetic = -hbar ** 2 / (2 * m) * d2psi
        potential = self.V * psi

        E_kin = _numeric_integral(np.conj(psi) * kinetic, self.x)
        E_pot = _numeric_integral(np.conj(psi) * potential, self.x)

        return np.real(E_kin + E_pot)

    def step(self) -> tuple[float, float]:
        """
        Advance the wave function by one time step.

        Returns:
            Tuple of (norm, energy_expectation) after the step
        """
        # Crank-Nicolson: A @ psi^{n+1} = B @ psi^n
        rhs = self.B @ self.psi
        self.psi = self.spsolve(self.A, rhs)

        # Compute diagnostics
        norm = self.compute_norm(self.psi)
        energy = self.compute_energy_expectation(self.psi)

        return norm, energy

    def run(self, callback: Optional[Callable[[int, np.ndarray, float], None]] = None) -> SimulationResult:
        """
        Run the full time evolution.

        Args:
            callback: Optional function called each frame: callback(frame, psi, norm)

        Returns:
            SimulationResult containing all time steps
        """
        num_steps = int(self.config.total_time / self.config.time_step)

        # Storage arrays
        prob_density = np.zeros((num_steps, len(self.x)))
        real_part = np.zeros((num_steps, len(self.x)))
        imag_part = np.zeros((num_steps, len(self.x)))
        self.norm_history = np.zeros(num_steps)
        self.energy_history = np.zeros(num_steps)

        time_axis = np.linspace(0, self.config.total_time, num_steps)

        print(f"Running simulation: {num_steps} time steps...")

        for j in range(num_steps):
            # Store current state
            prob_density[j] = np.abs(self.psi) ** 2
            real_part[j] = np.real(self.psi)
            imag_part[j] = np.imag(self.psi)

            # Advance and record diagnostics
            norm, energy = self.step()
            self.norm_history[j] = norm
            self.energy_history[j] = energy

            # Progress callback
            if callback is not None and j % max(1, num_steps // 100) == 0:
                callback(j, self.psi, norm)

            # Warn if norm increases (CAP should only decrease norm via absorption)
            if norm > 1.01:
                warnings.warn(f"Norm increased at step {j}: {norm:.6f} (possible instability)")

        print(f"Simulation complete. Final norm: {self.norm_history[-1]:.6f}")
        print(f"Energy conservation: {self.energy_history[0]:.6f} -> {self.energy_history[-1]:.6f}")

        return SimulationResult(
            x=self.x,
            time_axis=time_axis,
            probability_density=prob_density,
            real_part=real_part,
            imag_part=imag_part,
            norm_history=self.norm_history,
            energy_history=self.energy_history
        )


# =============================================================================
# Visualization with Interactive Buttons
# =============================================================================

class InteractiveSimulation:
    """Interactive Schrödinger simulation with preset buttons."""

    # Preset configurations
    PRESETS = {
        "Harmonic Trap": {
            "config": SimulationConfig(
                x_min=0.0, x_max=10.0, num_points=300, total_time=8.0,
                time_step=0.001, packet_center=3.5, packet_width=0.4,
                packet_momentum=0.0, cap_enabled=True, cap_width=1.0, cap_strength=0.3,
            ),
            "potential": lambda x: harmonic_oscillator(x, k=20.0, center=5.0),
            "initial_state": lambda x, cfg: gaussian_wave_packet(x, cfg.packet_center, cfg.packet_width, cfg.packet_momentum),
            "description": "Particle oscillating in harmonic trap",
        },
        "Gaussian Well": {
            "config": SimulationConfig(
                x_min=0.0, x_max=10.0, num_points=300, total_time=8.0,
                time_step=0.001, packet_center=4.5, packet_width=0.4,
                packet_momentum=0.0, cap_enabled=True, cap_width=1.0, cap_strength=0.3,
            ),
            "potential": lambda x: gaussian_well(x, depth=100.0, center=5.0, width=2.0),
            "initial_state": lambda x, cfg: gaussian_wave_packet(x, cfg.packet_center, cfg.packet_width, cfg.packet_momentum),
            "description": "Particle oscillating in deep Gaussian well",
        },
        "Finite Well": {
            "config": SimulationConfig(
                x_min=0.0, x_max=10.0, num_points=300, total_time=8.0,
                time_step=0.001, packet_center=5.0, packet_width=0.4,
                packet_momentum=-15.0, cap_enabled=True, cap_width=1.0, cap_strength=0.4,
            ),
            "potential": lambda x: np.where((x > 3.0) & (x < 7.0), -500.0, 0.0),
            "initial_state": lambda x, cfg: gaussian_wave_packet(x, cfg.packet_center, cfg.packet_width, cfg.packet_momentum),
            "description": "Particle trapped in finite square well",
        },
    }

    def __init__(self):
        self.fig = None
        self.ax = None
        self.current_animation = None

    def run_preset(self, preset_name: str):
        """Run a specific preset simulation."""
        preset = self.PRESETS[preset_name]
        config = preset["config"]
        potential_func = preset["potential"]
        initial_state_func = preset["initial_state"]

        print(f"\nRunning: {preset_name}")
        print(f"Description: {preset['description']}")

        # Create solver and run
        initial_state = lambda x: initial_state_func(x, config)
        solver = SchrodingerSolver(config, potential_func, initial_state)
        result = solver.run()

        # Update or create plot
        if self.fig is None:
            self._create_plot()

        self._update_plot(result, config, potential_func, preset_name)

    def _create_plot(self):
        """Create the figure and axes."""
        self.fig, self.ax = plt.subplots(figsize=(10, 7))
        plt.subplots_adjust(bottom=0.25)

        # Create button axes
        button_y = 0.02
        button_height = 0.04
        n_buttons = len(self.PRESETS)
        button_width = 0.18
        spacing = 0.02

        # Calculate starting x to center buttons
        total_width = n_buttons * button_width + (n_buttons - 1) * spacing
        start_x = (1 - total_width) / 2

        self.buttons = []
        for i, name in enumerate(self.PRESETS.keys()):
            ax_btn = plt.axes([start_x + i * (button_width + spacing), button_y, button_width, button_height])
            btn = Button(ax_btn, name, color='lightgoldenrodyellow', hovercolor='yellow')
            btn.on_clicked(lambda event, n=name: self._run_preset_and_close(n))
            self.buttons.append(btn)

        self.ax.set_xlabel("x")
        self.ax.set_ylabel("|ψ|²")
        self.ax.grid(True, alpha=0.3)

    def _run_preset_and_close(self, preset_name: str):
        """Handle button click by closing current figure and running preset."""
        plt.close(self.fig)
        self.fig = None
        self.run_preset(preset_name)

    def _update_plot(self, result: SimulationResult, config: SimulationConfig,
                     potential_func: Callable, preset_name: str):
        """Update the plot with new simulation data."""
        self.ax.clear()
        self.ax.set_xlim(config.x_min, config.x_max)
        self.ax.set_ylim(0, max(0.3, np.max(result.probability_density) * 1.5))
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("|ψ|²")
        self.ax.set_title(f"Code of Zain Ul Abideen\nUpgraded with Claude Code + Qwen3.5\n\n{preset_name}")
        self.ax.grid(True, alpha=0.3)

        # Plot potential (scaled to fit)
        V = potential_func(result.x)
        V_min, V_max = np.min(V), np.max(V)
        if V_max > V_min:
            V_scaled = (V - V_min) / (V_max - V_min + 1e-10) * (self.ax.get_ylim()[1] * 0.8)
        else:
            V_scaled = V
        self.line_pot, = self.ax.plot(result.x, V_scaled, 'r-', linewidth=2, alpha=0.8, label="Potential")

        # Plot probability density
        self.line_prob, = self.ax.plot([], [], 'b-', linewidth=2, label="|ψ|²")
        self.ax.legend(loc='upper right')

        # Info text
        self.info_text = self.ax.text(0.02, 0.98, '', transform=self.ax.transAxes,
                                       fontsize=10, verticalalignment='top',
                                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        def init():
            self.line_prob.set_data([], [])
            self.info_text.set_text('')
            return self.line_prob, self.info_text

        def update(frame):
            self.line_prob.set_data(result.x, result.probability_density[frame])
            self.info_text.set_text(f'Time: {result.time_axis[frame]:.3f}')
            return self.line_prob, self.info_text

        self.current_animation = FuncAnimation(self.fig, update, frames=len(result.time_axis),
                                                init_func=init, blit=True, interval=10)
        plt.show()


def main():
    """Main entry point with interactive buttons."""
    app = InteractiveSimulation()

    # Start with harmonic trap by default
    app.run_preset("Harmonic Trap")


if __name__ == "__main__":
    main()
