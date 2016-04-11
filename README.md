# phd-thesis
Source code and data for the thesis "Accelerating Monte Carlo methods for Bayesian inference in dynamical models"

The material will be uploaded before May 4.

## Example: Reconstructing the temperature of pre-historic Earth
The code in the folder **example-icevarve** replicates the example in Section 1.1.1. Two different models are used: (i) a Gaussian process regression and (ii) a non-linear state space model SSM.

For (i), we make use of a zero mean function and a kernel consisting of a bias kernel and a Matérn 5/2 kernel. The hyperparameters in the prior are hand-tuned and based on a run of empirical Bayes. The computations are carried out by **ex-icevarve-gp.py** and the plotting is done using **ex-icevarve-plot.R**. 

For (ii), we make use of a state space model proposed by Langrock (2011) < http://www.tandfonline.com/doi/abs/10.1080/02664763.2011.573543 > given by x_{t+1} \sim \mathcal{N}(x_{t+1}; \phi x_t, \sigma^2_v), y_t \sim \mathcal{G}(y_t; \alpha, \beta \exp(-x_t)), where \{\phi,\sigma_v,\alpha,\beta\} denote the parameters. These are estimated using the particle Metropolis-Hastings algorihtm by running the file **ex-icevarve-pmh0.py**, which also estimates the latent state by marginalising over the parameter posterior. The results are plotted using **ex-icevarve-plot.R**.

## Example: How does unemployment affect inflation?

## Example: Voting behaviour in the US Supreme court
