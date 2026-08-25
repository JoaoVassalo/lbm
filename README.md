[🇧🇷 Português](README.pt-br.md)


# LBM-CUDA

Two-dimensional Lattice Boltzmann Method (LBM) solver developed in CUDA, using moment-based storage and parallel execution on the GPU.

## Structure

- `src/config/` — geometry, D2Q9 stencil, moments, physics, and time parameters.
- `src/lbm/init/` — domain and lattice initialization.
- `src/lbm/simulation/` — streaming, distribution reconstruction, and collision.
- `src/io/` — VTK output and performance measurements.
- `src/core/` — indexing and auxiliary types.
- `run.sh` — compiles and runs the solver.

## Method

The code uses the **D2Q9** stencil and stores six variables per node:

\[
(\rho, u_x, u_y, m_{xx}, m_{xy}, m_{yy}).
\]

The distributions are reconstructed from these moments, followed by the streaming and collision steps. The simulation uses two buffers (`momA`/`momB`) in a ping-pong scheme.

## Execution

`nvcc` must be available in the `PATH`.

```bash
bash run.sh