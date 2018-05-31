def Gaussmixprob(x, n_Gauss, mus, sigmas, amps):
    # Computes the probability prob of an event x if it is part of a Gaussian mixture distribution with n_Gauss number of components
    import math
    component_sum = 0
    normalisation_sum = 0
    for i in range(n_Gauss):
        component = amps[i]*math.exp(-((x-mus[i])**2/(2*(sigmas[i]**2))))
        component_sum += component 
        component_norm_factor = amps[i]*sigmas[i]
        normalisation_sum += component_norm_factor
    prob = ((math.sqrt(2*math.pi)*normalisation_sum)**-1)*component_sum
    return prob
    
def fit(data, n_Gauss, ini_mu, ini_sigma, ini_amp, tol, max_it):
    # ===========================================================================
    # Uses an expectation-maximisation algorithm to fit a Gaussian mixture model
    # with "n_Gauss" components to a 1D data set. Include initial guesses (ini_mu,
    # ini_sigma, ini_amp) for the fitted parameters and a tolerance value for the
    # change in log-likelihood at which the algorithm stops iterating. Max_it
    # specifies the maximum number of iterations performed.
    # ===========================================================================
    import numpy as np
    import math
    
    def normprob(x, mu, sigma):
        # Computes the probability prob of an event x if it's part of the normal distribution with mu and sigma
        prob = (1/(math.sqrt(2*math.pi*(sigma**2))))*math.exp(-((x-mu)**2/(2*(sigma**2))))
        return prob
    
    def bayesprob(x, case, n_Gauss, mus, sigmas, amps):
        # Computes the probability that x is part of a Gaussian distribution defined by the index in "case"
        # out of a total number of Gaussians (n_Gauss with mus, sigmas and relative proportions amps)
        prob = np.zeros(n_Gauss)
        for i in range(n_Gauss):
            prob[i] = normprob(x, mus[i], sigmas[i])*amps[i]

        posterior = prob[case]/prob.sum()
        return posterior
    
    def mu_estimator(data, posteriors):
        # Estimates a mean for the data given certain posterior probabilities
        return (posteriors*data).sum()/posteriors.sum()
    
    def sigma_estimator(data, posteriors, mu):
        # Estimates a sigma for the data given certain posterior probabilities and mean (mu)
        weighted_variances = posteriors*(data-mu)**2
        variance_estimate = weighted_variances.sum()/posteriors.sum()
        sigma_estimate = math.sqrt(variance_estimate)
        return sigma_estimate
    
    def amp_estimator(posteriors):
        # Estimates a relative amplitude for a Gaussian given posterior probabilities for a dataset
        amp_estimate = posteriors.sum()/float(len(posteriors))
        return amp_estimate
    
    
    # Start of actual algorithm
    data = np.asarray(data)
    mus = np.asarray(ini_mu)
    sigmas = np.asarray(ini_sigma)
    amps = np.asarray(ini_amp)

    # Compute log-likelihood at start of optimisation
    likelihoods = np.zeros(len(data))
    for i in range(len(data)):
        likelihoods[i] = (Gaussmixprob(data[i], n_Gauss, mus, sigmas, amps))
    log_likelihoods = np.log(likelihoods)

    log_likelihood_start = np.sum(log_likelihoods)
    n = 0
    while n <= max_it:
        # E-step
        # Calculate posterior probability for all data points assuming they are part of one of the Gaussians

        posteriors = np.zeros((n_Gauss,len(data)))
        for i in range(n_Gauss):
            for j, point in enumerate(data):
                posteriors[i,j] = bayesprob(point, i, n_Gauss, mus, sigmas, amps)
           
        # M-step
        # Calculate better estimates for mus, sigmas and proportions based on the posterior probabilities

        mus = np.zeros(n_Gauss)
        sigmas = np.zeros(n_Gauss)
        amps = np.zeros(n_Gauss)

        for i in range(n_Gauss):
            mus[i] = mu_estimator(data, posteriors[i])
            sigmas[i] = (sigma_estimator(data, posteriors[i], mus[i]))
            amps[i] = amp_estimator(posteriors[i])

        # Compute log-likelihood at end of iteration step and how it's changed during the iteration step
        likelihoods = np.zeros(len(data))
        for i in range(len(data)):
            likelihoods[i] = (Gaussmixprob(data[i], n_Gauss, mus, sigmas, amps))
        log_likelihoods = np.log(likelihoods)

        log_likelihood = np.sum(log_likelihoods)

        delta_loglike = abs(log_likelihood - log_likelihood_start)
    
        log_likelihood_start = log_likelihood
        
        if n == max_it:
            print('Warning, maximum iterations reached')
        
        n += 1
        
        if delta_loglike > tol:
            continue
        else:
            break
    
    return mus, sigmas, amps, n, log_likelihood