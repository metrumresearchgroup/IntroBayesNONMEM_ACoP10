# IntroBayesNONMEM

## Introduction to Bayesian pharmacometric data analysis using NONMEM

This repo contains the developing content for an introductory workshop on Bayesian pharmacometric data analysis using NONMEM.

### Objectives
* Introduce principles and methods of Bayesian data analysis for pharmacometric applications.
* Provide a hands-on experience in Bayesian data analysis using NONMEM.

### Primary intended audience: Pharmacometricians
### Background assumed
* Population PKPD modeling
* Use of R and NONMEM

## Workshop outline
* Why Bayesian?
* Introduction to Bayesian statistical principles and methods
  * Bayes Rule
  * Bayesian modeling \& inference process
* Computation for Bayesian modeling
  * Maximum a Posteriori (MAP) Bayes
    * Individual: NONMEM POSTHOC
    * Population: Penalized Maximum Likelihood
  * Full Bayesian analysis
    * General computational approach: posterior simulation
    * Brief intro to Markov chain Monte Carlo (MCMC) simulation
      * Gibbs sampling 
      * Metropolis-Hastings 
      * Hamiltonian Monte Carlo and NUTS
* Overview of NONMEM implementations
  * MAP estimation
    * Using prior distributions with optimization methods
  * MCMC: BAYES and NUTS methods
  * Prior specification in NONMEM
* Hands-on 1: Example illustrating Bayesian data analysis workflow
* Prior distributions
  * Role of a prior distribution
  * Informative, uninformative or weakly informative?
* Hands-on 2: MAP popPK with selective use of informative priors for nuisance parameters: pediatric atorvastatin
* Model evaluation and comparison
* Assessing convergence and choosing numbers of burn-in and post-burn-in samples
* Getting your hands on posterior samples for individual parameters and predictions
* Hands-on 3: Full Bayes popPK with selective use of informative priors for nuisance parameters: pediatric atorvastatin
* When stuff goes wrong
  * Diagnosing and remedying sampling problems encountered with MCMC
  * Reparameterization, e.g., centered vs non-centered parameterizations for hierarchical models
  * Prior distributions as part of the solution
* Hands-on 4: Full Bayes popPKPD using semi-mechanistic model
  * Friberg-Karlsson semi-mechanistic model for drug-induced myelosuppression
  * Informative priors for drug-independent system parameters
* Practical strategies for selecting Bayesian estimation methods for specific types of problems
  * When to go Bayes (and why)?
  * Which method?
  * Which tool?
* Preview of Bayesian data analysis using Stan and Torsten
  * Brief intro with demo
  * Advantages/disadvantages
* What didn't we cover?
