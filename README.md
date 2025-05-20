**Bayesian Random Phase-Amplitude Gaussian Process (RPAGP)**

**Overview**

This project implements the Random Phase-Amplitude Gaussian Process (RPAGP) model, a flexible Bayesian approach to estimating trial-level Event-Related Potentials (ERPs) from EEG data. Unlike traditional ERP methods that average over trials and obscure individual differences, RPAGP captures trial-specific amplitude and latency variations while modeling background brain activity with an autoregressive process.

The project also addresses the computational challenges of Gaussian Process (GP) models by incorporating two acceleration techniques: the projection method and the nearest-neighbor method.

⸻

**Objectives**
	•	Improve single-trial ERP estimation in EEG data using Bayesian modeling
	•	Implement the RPAGP model with full posterior inference via MCMC
	•	Reduce computational cost with scalable GP approximations
	•	Evaluate performance on both empirical and simulated datasets

⸻

**Key Features**
	•	Flexible GP modeling of the ERP waveform
	•	Bayesian inference for trial-level amplitude and latency shifts
	•	Autoregressive noise modeling for structured background EEG activity
	•	Projection method and Nearest-Neighbor GP for computational efficiency
	•	Signal-to-noise ratio (SNR) and standardized measurement error (SME) metrics for model evaluation
	•	Comparison of empirical vs. model-based ERP estimation

⸻

**Dataset**
	•	161 participants from a cognitive neuroscience experiment
	•	ERP response to emotionally charged images (high, low, and neutral arousal)
	•	EEG recordings segmented into 900-ms epochs using 129-channel sensor net

⸻

**Modeling Highlights**
	•	Observed ERP signal for each trial is modeled as a transformation of a shared latent signal:
y_i(t) = β_i * f(t - τ_i) + v_i(t)
	•	Shared signal f modeled via Gaussian Process
	•	Trial-specific noise v_i(t) follows an autoregressive (AR) process

⸻

**Computational Improvements**

**Projection Method**
	•	Projects full GP onto a smaller set of “knots”
	•	Reduces time complexity from O(n³) to O(m³)
	•	Achieves ~30–60x speedup with minimal loss in accuracy (evaluated using MSE)

**Nearest-Neighbor GP**
	•	Approximates GP by conditioning only on local neighbors
	•	Produces sparse precision matrices and allows parallel inference
	•	Scalable to massive spatial or time-series data

⸻

**Implementation**
	•	Language: Python
	•	Libraries: PyMC3, NumPy, SciPy, Matplotlib
	•	Inference via Metropolis-within-Gibbs MCMC
	•	Sensitivity tests confirm robustness to hyperparameter settings
	•	Run time: ~1–2 hours per subject (depending on size)

⸻

**Results**
	•	RPAGP outperforms traditional averaging in signal clarity, uncertainty quantification, and condition discrimination
	•	Model sharply narrows confidence intervals and improves detection of condition-level ERP differences
	•	Simulation shows projection methods retain high accuracy with dramatically lower computation time 
