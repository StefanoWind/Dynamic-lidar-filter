# Dynamic-lidar-filter
Dynamic filter of lidar data based on joint pdf of normalized radial wind speed and SNR (based on Beck and Kuhn, Remote sensing 9(6),  2017). 

# Inputs (lidar data)
 - t: timestamp, unit = epoch time, dimension = $N_t \times 1$
 - r: range, unit = m, dimension = $N_r \times 1$
 - ele: elevation, unit = $^\circ$, dimension   $N_t \times 1$
 - azi: azimuth, unit = $^\circ$, dimension   $N_t \times 1$
 - rws: radial wind speed, unit = m/s, dimension = $N_r \times N_t$
 - snr: signal-to-noise ratio, unit = dB, dimension = $N_r \times N_t$
 
# Inputs (filter)
 - Dx, Dy, Dz: bin size in $x,y,z$, unit = m, dimension = $1 \times 1$
 - Dt: bin size in time, unit = s, dimension = $1 \times 1$ (600 s recommended)
 - max_rws: max aboslute value of $u_{LOS}$, unit = m/s, dimension = $1 \times 1$ (30 m/s recommended)
 - min_p: min probability (normalized by the maximum), unit = N/A, dimension = $1 \times 1$ (0.0025 recommended)
 - min_N: min number of points per bin, unit = N/A, dimension = $1 \times 1$ (10 recommended)
 - min_ratio: min ratio of good points per bin, unit = N/A, dimension = $1 \times 1$ (0.2 recommended)

# Outputs (basic)
- rws_qc1: radial velocity after application of max_rws and min_p, unit = m/s, dimension = $N_r \times N_t$
- rws_qc2: radial velocity after application of all filters, unit = m/s, dimension = $N_r \times N_t$

# Outputs (diagnostic)
- rws_norm: normalized radial wind speed, unit = m/s, dimension = $N_r \times N_t$
- snr_norm: normalized signal-to-noise ratio, unit = dB, dimension = $N_r \times N_t$
- p_all: probability of $u^\prime_{LOS} - SNR^\prime$ (normalized by the maximum) for each data point, unit = N/A, dimension = $N_r \times N_t$
- N_all: number of data point per bin in $x,y,z,t$ for each data point, unit = N/A, dimension = $N_r \times N_t$
- ratio_all: ratio of good data point per bin in $x,y,z,t$ for each data point, unit = N/A, dimension = $N_r \times N_t$
- t_calc: computational time, unit = s, dimension = $12 \times 1$
- steps: computational steps, unit = string, dimension = $12 \times 1$

