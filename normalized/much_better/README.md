# Gaussian potential variant

This version uses an inverted Gaussian potential with a Gaussian initial wave
packet instead of the quadratic potential used in the parent folder. The wave
packet survives much longer, but the explicit finite-difference scheme still
becomes unstable after a long run.

The Crank–Nicolson solver in [`../../crank-nicolson/`](../../crank-nicolson/)
avoids this instability.
