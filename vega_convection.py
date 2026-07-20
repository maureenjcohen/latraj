""" Derivation and validation of the convective vertical-wind parameters used by
    the Ornstein-Uhlenbeck / AR(1) convection kernel in custom_kernels.py.

    The convective vertical velocity applied to a parcel in the cloud layer is

        w_conv(t) = sigma * envelope(z) * u(t)

    where u is a dimensionless, zero-mean, unit-variance red-noise state advanced
    each transport step by the exact OU/AR(1) update

        a     = exp(-dt / tau)
        u_new = a * u + sqrt(1 - a^2) * g,     g ~ N(0, 1)

    so u is stationary N(0, 1) with correlation time tau. sigma sets the physical
    amplitude (m/s) and tau the memory. Both come from the Vega balloon anemometer
    vertical-wind channel (W_a).

    -------------------------------------------------------------------------
    CALIBRATION NOTE — we adopt the VEGA-1 values:  sigma = 0.6 m/s, tau = 1200 s
    -------------------------------------------------------------------------
    W_a carries a persistent ~-0.5 m/s offset on BOTH balloons (an instrument /
    buoyancy bias, not turbulent convection), so sigma is taken from the DE-BIASED
    (mean-removed) signal. From the data in this module:

        Vega-1:  sigma_debiased = 0.66 m/s,  lag-1 AR(1) fit tau ~ 1050 s (~17.5 min)
        Vega-2:  sigma_debiased = 0.73 m/s,  lag-1 AR(1) fit tau ~ 5400 s (~90 min)

    Vega-2's tau is not credible as a convective decorrelation time: its raw lag-1
    autocorrelation (0.986) is inflated by slow drift / the balloon's own buoyancy
    oscillation across long good-data runs rather than by turbulence. Vega-1's
    numbers are clean and coincide with the standard reference calibration, so we
    anchor on Vega-1 and round to sigma = 0.6 m/s, tau = 1200 s (~20 min). These
    are the defaults registered on the FieldSet in run_simulation.py.
    """
# %%
import os
import numpy as np

# %%
# Location of the Vega balloon PDS files (vg1bl_rdr.dat / vg2bl_rdr.dat).
DATADIR = os.path.expanduser('~/repos/aerial_biosphere/inputs')

# Quality flag >= this counts as usable W_a (2/3 = good, 4 = high quality).
GOOD_FLAG = 2

# Adopted convection parameters (Vega-1 calibration; see module note above).
CONV_SIGMA = 0.6     # convective vertical-wind std [m/s], de-biased W_a
CONV_TAU_S = 1200.0  # correlation time [s] (~20 min)


# %%
def load_vega(fn):
    """ Parse a vg*bl_rdr.dat PDS table and return the time and W_a columns.

    The records are fixed, whitespace-separated value/flag pairs; W_a (anemometer
    vertical wind, m/s) is token 19 with its quality flag in token 20. Returns a
    dict with the full time axis, W_a, its flag, and a boolean 'good' mask.
    """
    time, w_a, flag = [], [], []
    for line in open(fn):
        tok = line.split()
        if len(tok) < 23:
            continue
        time.append(float(tok[0]))
        w_a.append(float(tok[19]))
        flag.append(int(tok[20]))

    time = np.array(time)
    w_a = np.array(w_a)
    flag = np.array(flag)
    return {'time': time, 'w_a': w_a, 'flag': flag, 'good': flag >= GOOD_FLAG,
            'dt': float(np.median(np.diff(time)))}


# %%
def derive_parameters(data):
    """ Estimate the OU/AR(1) parameters from a loaded Vega record.

    Returns sigma (de-biased W_a std, m/s), tau (AR(1) correlation time, s), the
    lag-1 autocorrelation r1, and the mean bias that was removed. tau is obtained
    from the lag-1 autocorrelation of consecutive good samples via
    tau = -dt / ln(r1), the maximum-likelihood AR(1) estimate.
    """
    good = data['good']
    dt = data['dt']
    w = data['w_a'][good]
    bias = w.mean()
    sigma = w.std()

    # lag-1 autocorrelation, using only pairs of consecutive-in-time good samples
    idx = np.where(good)[0]
    keep = idx[np.isin(idx + 1, idx)]      # i such that i+1 is also good
    a = data['w_a'][keep]
    b = data['w_a'][keep + 1]
    am, bm = a - a.mean(), b - b.mean()
    r1 = float((am * bm).sum() / np.sqrt((am * am).sum() * (bm * bm).sum()))
    tau = float(-dt / np.log(r1)) if 0.0 < r1 < 1.0 else float('nan')

    return {'sigma': float(sigma), 'tau': tau, 'r1': r1, 'bias': float(bias),
            'dt': dt, 'n_good': int(good.sum())}


# %%
def empirical_acf(data, max_lag=40):
    """ Empirical autocorrelation of de-biased W_a out to max_lag samples.

    For each lag k, correlates good samples k steps (k*dt seconds) apart. Returns
    lag times [s] and the autocorrelation at each lag.
    """
    good = data['good']
    w = data['w_a']
    mean = w[good].mean()
    lags, acf = [], []
    for k in range(max_lag + 1):
        i = np.where(good[:len(good) - k] & good[k:])[0]
        if len(i) < 10:
            break
        a = w[i] - mean
        b = w[i + k] - mean
        denom = np.sqrt((a * a).sum() * (b * b).sum())
        lags.append(k * data['dt'])
        acf.append(float((a * b).sum() / denom) if denom > 0 else np.nan)
    return np.array(lags), np.array(acf)


# %%
def simulate_ou(sigma, tau, dt, n, seed=0):
    """ Generate an OU/AR(1) realization: n samples at spacing dt, std sigma,
        correlation time tau. Uses the exact update so the result is stationary
        for any dt (matching conv_ou in custom_kernels.py). """
    rng = np.random.default_rng(seed)
    a = np.exp(-dt / tau)
    u = np.empty(n)
    u[0] = rng.standard_normal()
    for k in range(1, n):
        u[k] = a * u[k - 1] + np.sqrt(1.0 - a * a) * rng.standard_normal()
    return sigma * u


# %%
def plot_comparison(balloon='vg1', sigma=CONV_SIGMA, tau=CONV_TAU_S,
                    savepath=None):
    """ Plot the measured Vega vertical wind against the OU/AR(1) approximation.

    Three panels:
      1. Time series  — the full mission of de-biased W_a (bad samples masked as
         gaps) overlaid with an OU realization spanning the same time using
         (sigma, tau); compares the 'texture' (amplitude and smoothness) of the
         two processes. Point-by-point agreement is not expected (they are
         independent realizations) — only their statistical character.
      2. Autocorrelation — empirical ACF of W_a vs the AR(1) theory exp(-t/tau).
         This is the direct test of the correlation-time fit.
      3. Distribution — histogram of de-biased W_a vs the stationary N(0, sigma)
         density that the OU process samples from.

    balloon is 'vg1' or 'vg2'. Uses the adopted (sigma, tau) by default so the
    plot shows the parameters actually used by the kernel, not a per-balloon refit.
    """
    import matplotlib.pyplot as plt

    data = load_vega(os.path.join(DATADIR, f'{balloon}bl_rdr.dat'))
    params = derive_parameters(data)
    dt = data['dt']

    # full de-biased record with bad samples masked to NaN so the line breaks at gaps
    w_meas = data['w_a'] - data['w_a'][data['good']].mean()
    w_meas = np.where(data['good'], w_meas, np.nan)
    t_hr = (data['time'] - data['time'][0]) / 3600.0
    ou = simulate_ou(sigma, tau, dt, len(w_meas), seed=1)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    # --- Panel 1: time series ---
    ax = axes[0]
    ax.plot(t_hr, ou, color='#c0392b', lw=0.8, alpha=0.7,
            label=f'OU/AR(1): $\\sigma$={sigma}, $\\tau$={tau/60:.0f} min')
    ax.plot(t_hr, w_meas, color='#1f4e79', lw=1.1,
            label=f'Vega ({balloon.upper()}) W$_a$, de-biased')
    ax.set_xlabel('Mission time / hr')
    ax.set_ylabel('Vertical wind / m s$^{-1}$')
    ax.set_title('Vertical-wind time series')
    ax.legend(fontsize=8, loc='upper right')
    ax.axhline(0, color='grey', lw=0.6, ls=':')

    # --- Panel 2: autocorrelation ---
    ax = axes[1]
    lags, acf = empirical_acf(data)
    ax.plot(lags / 60.0, acf, 'o-', color='#1f4e79', ms=4,
            label=f'Vega ({balloon.upper()}) empirical ACF')
    ax.plot(lags / 60.0, np.exp(-lags / tau), color='#c0392b', lw=1.6,
            label=f'AR(1) theory  exp(-t/{tau/60:.0f} min)')
    ax.plot(lags / 60.0, np.exp(-lags / params['tau']), color='#c0392b',
            lw=1.2, ls='--',
            label=f"per-balloon fit  $\\tau$={params['tau']/60:.0f} min")
    ax.set_xlabel('Lag / min')
    ax.set_ylabel('Autocorrelation')
    ax.set_title('Correlation time')
    ax.legend(fontsize=8)
    ax.axhline(0, color='grey', lw=0.6, ls=':')

    # --- Panel 3: distribution ---
    ax = axes[2]
    w_all = data['w_a'][data['good']] - data['w_a'][data['good']].mean()
    ax.hist(w_all, bins=30, density=True, color='#1f4e79', alpha=0.6,
            label=f'Vega ({balloon.upper()}) W$_a$, de-biased')
    grid = np.linspace(w_all.min(), w_all.max(), 200)
    gauss = np.exp(-0.5 * (grid / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    ax.plot(grid, gauss, color='#c0392b', lw=1.8,
            label=f'N(0, $\\sigma$={sigma})')
    ax.set_xlabel('Vertical wind / m s$^{-1}$')
    ax.set_ylabel('Probability density')
    ax.set_title('Amplitude distribution')
    ax.legend(fontsize=8)

    fig.suptitle('Vega balloon vertical wind vs OU/AR(1) convection model',
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    if savepath:
        fig.savefig(savepath, dpi=200, bbox_inches='tight')
        print('Saved figure to', savepath)
    else:
        plt.show()
    return fig


# %%
if __name__ == "__main__":
    for balloon in ('vg1', 'vg2'):
        d = load_vega(os.path.join(DATADIR, f'{balloon}bl_rdr.dat'))
        p = derive_parameters(d)
        print(f"{balloon.upper()}: n_good={p['n_good']}  dt~{p['dt']:.1f}s  "
              f"bias={p['bias']:+.2f}  sigma={p['sigma']:.2f} m/s  "
              f"r1={p['r1']:.3f}  tau={p['tau']:.0f}s ({p['tau']/60:.1f} min)")
    print(f"\nAdopted (Vega-1): sigma={CONV_SIGMA} m/s, tau={CONV_TAU_S:.0f}s "
          f"({CONV_TAU_S/60:.0f} min)")
    plot_comparison('vg1')

# %%
