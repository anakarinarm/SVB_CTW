# Gaussian wind sensitivity experiments

Rplus_G5d and Rminus_G5d copy Rplus and Rminus respectively, using the same
bathymetry, grid, hydrography, initial state, OBCS physics, mixing, executable,
60-s timestep, and 1200-s diagnostics. No production files were changed.

Wind envelope: exp[-(t-216000)^2/(2*72000^2)], nominal peak 0.1 N/m^2.
The original binary vector map (including buffer taper) is factored from its
largest nonzero record; only its scalar temporal multiplier changes.
The actual original sampled peak is 0.0893407270 N/m^2. The review includes
both this actual series and its rescaling to the common 0.1 N/m^2 peak.

Run duration: 30 days (2592000 s); frequency spacing: 0.033333 cpd.
Forcing record spacing stays 1800 s. MITgcm uses midpoint records and linear
interpolation, so the Gaussian's two samples bracketing day 2.5 give a tiny
plateau from day 2.489583 to 2.510417, peaking at 0.0999921878 N/m^2
(0.0078% below the analytic maximum). No cutoff is applied at day 2.5 or 5.
At day 5 the analytic stress is 0.0011108997 N/m^2 and continues decaying.

The 31-day cycle prevents end-of-run wraparound. Its last record is a startup
guard evaluated at -900 s, used to bracket t=0; the run stops before this guard
is reached again. All other records sample the full Gaussian at midpoint times.
This yields a continuous interpolated startup and avoids a half-amplitude
artifact at t=0. Floating-point underflow occurs only for negligible late tails.

Boundary T/S files repeat the exact production records to 1488 records, matching
the shared forcing cycle. Only their paths and record counts change, not values.
data.pkg, data.diagnostics and eedata are byte-identical copies of the templates.
configuration.diff shows every changed namelist line.

Spectra use the same 30-day, 60-s sampled interval, no taper, no mean subtraction,
no zero padding, and squared Fourier amplitudes with time in days. The original
comparison pulse is interpolated from the actual binary record amplitudes and
scaled to the common peak, then continued at zero without repeating the original
16-day cycle. The Gaussian spectrum uses MITgcm interpolation.
New/original power at 1 cpd: 0.00108883.
Power ratio summed over 0.8--1.2 cpd: 0.00109775.
The nonzero starting Gaussian (1.1% of peak) contributes a small high-frequency
tail in a finite run; it is not the spectrum of an infinite Gaussian.
The longer pulse also supplies a much larger total wind impulse: this tests
duration/shape at fixed peak, not fixed integrated stress.
The scientific comparison targets the dominant frequency of the bay-generated
response rather than its absolute amplitude. No impulse normalization is applied
to either simulation; the configured Gaussian and run inputs remain unchanged.

Time-integrated stress from t=0 to day 30, at an untapered active-domain point:

| Forcing | Impulse (N s/m^2) | Impulse (N day/m^2) |
| --- | ---: | ---: |
| Original binary, actual | 367.629076 | 0.004254966 |
| Original, common 0.1 peak | 411.491028 | 0.004762628 |
| 5-day Gaussian, applied | 18023.348490 | 0.208603570 |
| 5-day Gaussian, formula | 18023.360991 | 0.208603715 |

The applied Gaussian supplies 43.800101
times the impulse of the common-peak original pulse, or
49.025906 times that of the actual
production binary pulse. Applied integrals use exact trapezoidal integration
of all MITgcm interpolation knots and both endpoints. The formula integral uses
the Gaussian error function. The original pulse is continued at zero after its
decay, without periodic repetition. These are stress impulses per unit area,
not domain-integrated forces. In the buffers, multiply by the unchanged spatial
taper. wind_stress_impulse.csv also records signed zonal/meridional components
at the untapered point; verification.json stores the scalar integrals and ratios.

Plots: forcing_comparison.png/pdf. Data: forcing_timeseries.csv,
gaussian_key_times.csv, forcing_spectrum.csv, wind_stress_impulse.csv.
Checks: verification.json.

MPI decomposition requires exactly 16 ranks (4 x 4), within the 20-core limit.
The launcher runs only the two new cases and refuses directories with outputs.
Output is written to LaCie. Existing diagnostics imply about 2.58 TB per case
for 30 days, excluding dumps/checkpoints; allow approximately 5.3 TB for the pair.
LaCie has approximately 5.6 TB available at setup. A 40-day pair would need
6.88 TB for diagnostics alone, so 30 days was selected to fit the drive.
Runtime is approximately twice a 15-day case at the same throughput.

Launch after reviewing the figures:

```bash
bash CTW_SVB/scripts/run_west64_south72_g5d_30d_sequential.sh
```

No run has been launched by this setup.
