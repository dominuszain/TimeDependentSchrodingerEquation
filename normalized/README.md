# Normalised wave packet

This is the same explicit finite-difference calculation as `../main.sce`, with
the wave function renormalised after each step. Renormalisation delays the growth
of the solution but does not remove it, because the scheme is only conditionally
stable. The `much_better/` folder changes the potential as well.

For a stable solver, see [`../crank-nicolson/`](../crank-nicolson/).
