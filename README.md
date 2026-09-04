# latraj — Lagrangian trajectories in planetary atmospheres

This repository takes wind fields from the Venus Planetary Climate Model
(VPCM), reformats them for the [OceanParcels](https://oceanparcels.org)
Lagrangian particle framework, advects virtual particles through those winds, and
plots the resulting 3-D trajectories through the Venus cloud decks.

---

## The pipeline

The project runs in three stages, one script each:

| Stage | Script | What it does |
|-------|--------|--------------|
| 1. Preprocess | `create_nc.py` | Reads a raw VPCM output file and writes one Parcels-ready netCDF file per timestep (`U`, `V`, `W` on time/height/lat/lon). |
| 2. Simulate | `run_simulation.py` | Builds a Parcels `FieldSet` from the preprocessed files, seeds particles, and advects them, writing positions to a `.zarr`. |
| 3. Visualise | `visualisations.py` | Plots the output `.zarr` trajectories as 3-D Plotly figures. |

Two supporting modules:

- `custom_kernels.py` — the physics kernels Parcels applies to each particle each
  timestep. We add Smagorinsky diffusion, Venus-tuned convection, periodic boundary conditions, and a balloon motion kernel.
- `vega_convection.py` — derives and validates the Venus convection kernel's parameters
  from the Vega balloon data, and plots the fit (`plot_comparison`).

---

## Setup

1. **Create the Python environment.** The heavy dependency is OceanParcels
   (which also pulls in xarray, netCDF4, and a C compiler for its JIT kernels);
   plotting uses plotly, matplotlib, or pillow. Use the environment file, e.g.

   ```bash
   conda env create -f environment.yml
   conda activate latraj
   ```

2. **Create your local config.** Paths and per-run settings live in `config.py`,
   which is **git-ignored** so machine-specific values never end up in commits or
   pull requests. Copy the tracked template and edit it:

   ```bash
   cp config_example.py config.py
   # then edit config.py: set the input/output paths and, if you like,
   # the particle start positions and run length.
   ```

   `config_example.py` is the documented list of every setting the scripts
   expect — if a script complains `No module named 'config'`, you skipped this
   step.

---

## Running

```bash
python create_nc.py        # stage 1: preprocess VPCM output -> Parcels netCDFs
python run_simulation.py   # stage 2: run the trajectory simulation -> .zarr
```

For stage 3, open the output `.zarr` as an xarray dataset and call a plotting
function, e.g. in a notebook or an interactive session:

```python
import xarray as xr
from visualisations import traj3d, doubletraj
ds = xr.open_zarr('<SAVE_DIR>/<OUTPUT_NAME>')
traj3d(ds)          # all trajectories on one 3-D plot
```

The scripts use `# %%` cell markers, so they also run cell-by-cell in VS Code or
a Jupyter-style interactive window, which is handy while developing.

---

## Venus data conventions (read before writing kernels)

These are easy to get wrong and are baked into the preprocessing:

- **Venus rotates retrograde.** The raw VPCM longitude axis runs +180 to -180;
  `create_nc.py` reverses it to -180..+180, flips all data arrays to match, and
  negates the U-wind. Keeping the raw orientation caused immediate out-of-bounds
  errors in Parcels. See the docstring at the top of `create_nc.py`.
- **The vertical coordinate is altitude, and up is positive.** Unlike the usual
  OceanParcels "depth increases downward" convention, here increasing
  `depth` / `z` means higher altitude (metres). A positive vertical velocity
  moves a particle **up**. Add upward terms as positive `particle_ddepth`.
- **The winds are treated as an A-grid** (U/V/W co-located), even though the PCM
  runs on an Arakawa C grid — the outputs are regridded before writing.

---

## Writing kernels (Parcels gotchas)

Parcels compiles kernels to C ("JIT"), so kernel functions are more restricted
than normal Python. When adding to `custom_kernels.py`:

- A kernel has the signature `def name(particle, fieldset, time):` and returns
  nothing; it changes the particle by adding to the special accumulators
  `particle_ddepth`, `particle_dlat`, `particle_dlon` (displacements this step),
  or by assigning to particle variables.
- You **cannot call your own Python helper functions** from inside a kernel —
  inline the maths instead (see how the convection envelope is inlined in
  `convection_ou`).
- Use `math.` functions, not `numpy`; draw random numbers with
  `parcels.ParcelsRandom.normalvariate(...)`.
- Per-particle state needs a `Variable` on the particle class (see `u_conv` on
  `VenusParticle`); run-wide constants go on the `fieldset` via
  `fieldset.add_constant(...)`.

`convection_ou` and `smagdiff` in `custom_kernels.py` are worked examples to copy
from.

---

## Working with git

`main` is the protected, reviewed branch. Do your work on a branch and open a
pull request for review — don't commit to `main` directly.

```bash
git checkout main && git pull        # start from the latest main
git checkout -b my-feature           # make a branch for your work
# ... edit, then ...
git add <the files you changed>
git commit -m "Short description of what changed"
git push -u origin my-feature        # push your branch
# then open a Pull Request on GitHub for review
```

Guidelines:

- **Keep run-config changes out of your PRs.** Edit particle positions or
  paths in `config.py` (git-ignored), not in the scripts. Your PRs should
  be about new physics or plots.
- **Test before you submit** (see any tests provided, or add one alongside your
  change).
- Don't commit generated files — outputs (`.zarr`, `.nc`), figures (`.png`), and
  `__pycache__` are all git-ignored on purpose.
- One focused change per branch/PR makes review much easier.
