# WindDancer: Understanding Acoustic Sensing under Ambient Airflow

## About

This repository contains the official implementation of "WindDancer: Understanding Acoustic Sensing under Ambient Airflow" by Kuang Yuan et al.

**Abstract:** Acoustic sensing has recently garnered significant interest for a wide range of applications ranging from motion tracking to health monitoring. However, prior works overlooked an important real-world factor affecting acoustic sensing systems—the instability of the propagation medium due to ambient airflow. Airflow introduces rapid and random fluctuations in the speed of sound, leading to performance degradation in acoustic sensing tasks. This paper presents WindDancer, the first comprehensive framework to understand how ambient airflow influences existing acoustic sensing systems, as well as provides solutions to enhance systems performance in the presence of airflow. Specifically, our work includes a mechanistic understanding of airflow interference, modeling of sound speed variations, and analysis of how several key real-world factors interact with airflow. Furthermore, we provide practical recommendations and signal processing solutions to improve the resilience of acoustic sensing systems for real-world deployment. We envision that WindDancer establishes a theoretical foundation for understanding the impact of airflow on acoustic sensing, and advances the reliability of acoustic sensing technologies for broader adoption.

## Repository Overview

This repository provides:
- **Simulation Framework**: AR-modeled speed-of-sound variations for chirp-based FMCW acoustic sensing
- **Real-World Solutions**: Two novel signal processing methods demonstrated on real acoustic data
  - Spectral variance-based ranging for static clutter rejection
  - VMD-based micro-motion tracking for breathing detection
- **Complete Pipeline**: From raw measurements → AR model estimation → simulation/evaluation

## Repository Structure

```
.
├── Simulation Scripts
│   ├── soundspeed_ar_modeling.m            # Estimate AR parameters from measurements
│   ├── example_simple_simulation.m          # Educational demo (no external data required)
│   └── simuliation_chirp_parasbatch.m      # Main simulation with AR-modeled SoS variations
│
├── Real-World Evaluation Scripts
│   ├── ranging_solution_evaluation.m        # Spectral variance-based ranging solution
│   └── phase_tracking_solution_evaluation.m # VMD-based micro-motion tracking solution
│
├── Utility Functions (utils/)
│   ├── generate_ar_noise.m                 # AR noise sequence generation
│   ├── fine_grained_eval.m                 # Micro-motion accuracy evaluation
│   ├── generate_transmit_sw.m              # FMCW chirp signal generation
│   ├── mixing_sw.m                         # I/Q demodulation
│   ├── SpecVar_thresholding.m              # Spectral variance clutter rejection
│   └── VMD_complex.m                       # Variational mode decomposition
│
├── Data Directories
│   ├── soundspeed/                         # Speed-of-sound data
│   │   ├── ar_models/                      # Pre-computed AR models
│   │   │   ├── distances-50-230_ar_models.mat
│   │   │   └── chirp-lengths-30-200_ar_models.mat
│   │   └── raw_measurements/               # Raw speed-of-sound data
│   │       ├── distances-50-230.mat
│   │       ├── chirp-lengths-30-200.mat
│   │       └── nowind.mat
│   └── realworld_data/                     # Example audio recordings
│       ├── example_1.wav                   # For ranging evaluation
│       └── example_2.wav                   # For phase tracking evaluation
│
└── README.md
```

## Requirements

### Software
- MATLAB R2019b or later
- Signal Processing Toolbox (for `fir1`, `filtfilt`, `findpeaks`, `vmd`)

### Hardware
- Standard PC (simulation runs in 1-10 seconds per AR parameter set)
- Audio interface (for collecting new real-world data, optional)

## Quick Start

### 1. Run Simple Example (No External Data Required)
```matlab
run('example_simple_simulation.m')
```
This demonstrates basic FMCW acoustic sensing without requiring any data files.

**Optional: Generate Custom AR Models**
```matlab
% Edit soundspeed_ar_modeling.m to specify your measurement file
run('soundspeed_ar_modeling.m')
```
This is only needed if you have new speed-of-sound measurements to process.

### 2. Run Main Simulation (AR-Modeled Variations)
```matlab
run('simuliation_chirp_parasbatch.m')
```
Simulates ranging under speed-of-sound variations using pre-computed AR models.

### 3. Evaluate Real-World Ranging Performance
```matlab
run('ranging_solution_evaluation.m')
```
Demonstrates spectral variance-based clutter rejection on real acoustic data.

### 4. Evaluate Real-World Micro-Motion Tracking
```matlab
run('phase_tracking_solution_evaluation.m')
```
Demonstrates VMD-based breathing motion detection on real acoustic data.

## Key Contributions

### 1. AR Modeling of Speed-of-Sound Variations
- Models temporal correlations in airflow-induced speed-of-sound fluctuations
- Enables realistic simulation of environmental interference
- **Script**: `soundspeed_ar_modeling.m`, `simuliation_chirp_parasbatch.m`

### 2. Spectral Variance-Based Ranging
- Novel clutter rejection method for ranging in presence of static reflectors
- Distinguishes moving targets from static objects using spectral variance analysis
- **Script**: `ranging_solution_evaluation.m`

### 3. VMD-Based Micro-Motion Tracking
- Variational Mode Decomposition to isolate breathing signals from interference
- Robust phase-based displacement estimation
- **Script**: `phase_tracking_solution_evaluation.m`

## System Parameters

- **Sampling Frequency**: 48 kHz
- **Carrier Frequency**: 18 kHz (ultrasonic)
- **Bandwidth**: 4 kHz
- **Chirp Duration**: 50-150 ms (configurable)


## Data Files

The repository includes pre-computed AR models in `soundspeed/ar_models/`:
- `distances-50-230_ar_models.mat` - AR parameters for different distances from airflow source
- `chirp-lengths-30-200_ar_models.mat` - AR parameters for different chirp durations
- Each contains: `ar_paras` (coefficients + noise std), `dlist`/`clist` (conditions), `plist` (model orders)

Raw measurements in `soundspeed/raw_measurements/`:
- `distances-50-230.mat` - Speed-of-sound measurements at different distances
- `chirp-lengths-30-200.mat` - Measurements with different chirp durations
- `nowind.mat` - Baseline measurements without airflow

Real-world audio examples in `realworld_data/`:
- `example_1.wav` - Recording for ranging evaluation
- `example_2.wav` - Recording for breathing motion detection

## Generating New AR Models

1. Collect speed-of-sound measurements under your conditions
2. Save in `soundspeed/raw_measurements/` as MAT file with `sos_est_all` variable
3. Edit `soundspeed_ar_modeling.m` to set distance list and AR orders
4. Run script and uncomment save command

## Extending the Code

All scripts are extensively documented with:
- Comprehensive headers explaining methodology
- Inline comments for each processing step
- Usage examples in function documentation
- References to relevant literature

Modify parameters in clearly marked sections (e.g., "SYSTEM PARAMETERS", "SCENE CONFIGURATION").

## Citation

If you find our work useful in your research, please consider citing:

```bibtex
@article{yuan2025windDancer,
  title={WindDancer: Understanding Acoustic Sensing under Ambient Airflow},
  author={Yuan, Kuang and Li, Dong and Zhou, Hao and Li, Zhehao and Qiu, Lili and Kumar, Swarun and Xiong, Jie},
  journal={Proceedings of the ACM on Interactive, Mobile, Wearable and Ubiquitous Technologies},
  volume={9},
  number={2},
  year={2025},
  publisher={ACM New York, NY, USA}
}
```
