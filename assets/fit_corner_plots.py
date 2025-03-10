import numpy as np
import scipy.optimize
import corner
import emcee

import matplotlib.pyplot as plt

def plot_data(x, y, y_err, models=None):
    fig, ax = plt.subplots()

    ax.errorbar(x, y, y_err, fmt=".", label="Data")

    if models is not None:
        for model_def in models:
            if "lower" in model_def:
                ax.fill_between(
                    model_def["x"], model_def["lower"], model_def["upper"],
                    **model_def.get("style", {})
                )
            else:
                ax.plot(
                    model_def["x"], model_def["y"],
                    **model_def.get("style", {})
                )

    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.legend(frameon=False)

    return fig, ax


def analyse_data(data,
                 log_posterior_fn, model_fn, predict_fn,
                 param_names, theta_true):

    x, y, sigma_y = data["x"], data["y"], data["y_err"]

    theta_init = theta_true

    true_model_plot_def = dict(
        x=x, y=model_fn(*theta_true, x),
        style=dict(color="black", ls="--", label="True model")
    )

    # Find the MAP
    def negative_log_posterior(theta, x, sigma_y, y):
        return -log_posterior_fn(*theta, x, sigma_y, y)
    
    MAP_result = scipy.optimize.minimize(
        fun=negative_log_posterior,
        x0=theta_init,
        args=(x, sigma_y, y)
    )
    theta_MAP = MAP_result.x

    print("MAP results")
    for name, theta in zip(param_names, theta_MAP):
        print(f"{name}_MAP = {theta}")
    
    plot_data(
        x=x, y=y, y_err=sigma_y,
        models=[
            true_model_plot_def,
            dict(x=x, y=model_fn(*theta_MAP, x),
                 style=dict(color="C1", label="MAP model")),
        ]
    )

    # Sample posterior with emcee

    # emcee passes an array of values for the sampled parameters
    # This wrapper just splits the array theta into m and b
    def log_posterior_wrapper(theta, x, sigma_y, y):
        return log_posterior_fn(*theta, x, sigma_y, y)

    # emcee requires some extra settings to run
    n_param = len(theta_true) # Number of parameter we are sampling
    n_walker = 10             # Number of walkers. This just needs to be larger than 2*n_param + 1
    n_step = 5000             # How many steps each walker will take. The number of samples will be n_walker*n_step

    # The starting point for each walker
    theta_init = theta_init + 0.1*np.random.normal(size=(n_walker, n_param))

    sampler = emcee.EnsembleSampler(
        nwalkers=n_walker, ndim=n_param,
        log_prob_fn=log_posterior_wrapper,
        args=(x, sigma_y, y)
    )
    state = sampler.run_mcmc(theta_init, nsteps=n_step, progress=True)

    # The samples will be correlated, this checks how correlated they are
    # We will discuss this once we come to MCMC methods
    print("Auto-correlation time of chain:")
    for name, value in zip(param_names, sampler.get_autocorr_time()):
        print(f"{name} = {value:.1f}")

    max_autocorr_time = max(sampler.get_autocorr_time())

    # We need to discard the beginning of the chain (a few auto-correlation times)
    # to get rid of the initial conditions
    chain = sampler.get_chain(
        discard=int(5*max_autocorr_time),
        thin=int(max_autocorr_time/2),
        flat=True
    )

    # Make a corner plot
    fig = plt.figure()
    fig = corner.corner(
        chain,
        bins=40,
        labels=param_names,
        truths=theta_true,
        levels=1-np.exp(-0.5*np.array([1, 2])**2), # Credible contours corresponding to 1 and 2 sigma in 2D
        quantiles=[0.025, 0.16, 0.84, 0.975],
        fig=fig
    )

    print("Posterior results (mean±std)")
    for i, name in enumerate(param_names):
        print(f"{name} = {np.mean(chain[:,i]):.2f}±{np.std(chain[:,i]):.2f}")

    # Make predictive distributions
    # Choose a small subsample of the chain for plotting purposes
    chain_samples = chain[np.random.choice(chain.shape[0], size=200)]
    # Evaluate the model at the sample parameters
    model_predictive = np.array(
        [model_fn(*sample, x) for sample in chain_samples]
    )
    model_quantiles = np.quantile(
        model_predictive, q=[0.025, 0.16, 0.84, 0.975], axis=0
    )

    # Get samples from the posterior predictive distribution
    posterior_predictive = np.array(
        [predict_fn(*sample, x, sigma_y) for sample in chain_samples]
    )
    posterior_predictive_quantiles = np.quantile(
        posterior_predictive, q=[0.025, 0.16, 0.84, 0.975], axis=0
    )

    # Make plots of the data, best-fit model, and predictive distributions
    plot_data(
        x=x, y=y, y_err=sigma_y,
        models=[
            true_model_plot_def,
            dict(x=x, y=model_fn(*theta_MAP, x),
                 style=dict(color="C1", label="MAP model")),
            
            dict(x=x, lower=model_quantiles[0], upper=model_quantiles[-1],
                 style=dict(color="C1", alpha=0.5, label="Model predictions")),
            dict(x=x, lower=model_quantiles[1], upper=model_quantiles[-2],
                 style=dict(color="C1", alpha=0.5)),

            dict(x=x, lower=posterior_predictive_quantiles[0], upper=posterior_predictive_quantiles[-1],
                 style=dict(color="grey", alpha=0.5, label="Posterior predictions")),
            dict(x=x, lower=posterior_predictive_quantiles[1], upper=posterior_predictive_quantiles[-2],
                 style=dict(color="grey", alpha=0.5)),
        ]
    )