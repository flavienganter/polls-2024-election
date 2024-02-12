# A Poll Aggregator for the 2024 French European Election

_Adaptation of the [poll aggregator implemented for the 2022 French Presidential Election](https://github.com/flavienganter/polls-2022-election)._

Estimates of the voting intentions for the 2024 French European Election (medians of the posterior distributions, and 95%, 90%, 80%, and 50% high density credible intervals):

![](https://github.com/flavienganter/polls-2024-election/blob/main/PollsFrance2024_latest.png?raw=true)

Evolutions of voting intentions since September 2023 (medians of the posterior distributions, and 95% and 50% high density credible intervals):

![](https://github.com/flavienganter/polls-2024-election/blob/main/PollsFrance2024_evolution.png?raw=true)

## Model

I use data from all voting intention polls fielded mid-June 2023, based on the survey reports available on the [Commission des sondages website](https://www.commission-des-sondages.fr/notices/). I build on [Heidemanns, Gelman and Morris (2020)](https://hdsr.mitpress.mit.edu/pub/nw1dzd02/release/1) to build a poll aggregator that does just that—aggregating polls—with no prediction intention whatsoever. The model is estimated with Stan.

For each scenario $i$ (part of poll $p\ =\ p_{[i]}$) and each candidate $c$, $s_{ci}^{\*}$ is the (adjusted) share of respondents who indicated support for candidate $c$. $s_{ci}^{\*}$ is typically not available, as polling firm round their estimates, so that one can only observe $s_{ci}$. To account for the uncertainty induced by the rounding, I model $s_{ci}^{\*}$ as a latent parameter defined by

$$ s_{ci}^{\*}\ =\ s_{ci}\ +\ \varepsilon_{ci} $$

with $\varepsilon_{ci}\ \sim\ \mathcal{U}\[b_{ci}^l;b_{ci}^u\]$, where $b_{ci}^l$ and $b_{ci}^u$ define the interval around $s_{ci}$ in which $s_{ci}^{\*}$ can be.

Noting $N_i$ the total number of respondents disclosing their voting intentions in the scenario $i$, I model the latent variable $s_{ci}^{\*}$ with a Beta distribution:

$$ s_{ci}^{\*}\ \sim\ \text{Beta}(\theta_{ci}N_i,(1-\theta_{ic})N_i) $$

where $\theta_{ci}$ is defined as

$$ \theta_{ci}\ \equiv\ \text{logit}^{-1}(\psi_c(date_i)\ +\ \mu_{cp\[i\]}\ +\ \lambda_{ch\[i\]}\  +\ X_i\beta_c(date_i)) $$

and $date_i$ is the date, centered so that $date_i\ =\ 1$ on June 16, 2023 (date of the beginning of the field period of the first survey considered).

### Splines

I model the evolutions of voting intentions over time with a spline of degree 3 with $K\ =\ 3$ knots:

$$ \psi_c(date_i)\ =\ \alpha_{c0}\cdot date_i\ +\ \sum_{k=1}^K\alpha_{ck}B_{3k}(date_i) $$

where $(B_{3k}(\cdot))\_k$ is a sequence of $B$-splines. I define a poll's date as the median day of the fielding period, or as the day immediately following the median when that median is not properly defined. To enforce smoothness and prevent the model from overfitting, I impose a random-walk prior on $(\alpha_{ck})\_{ck}$: $\alpha_{c0}\ \sim\ \mathcal{N}(0,\sigma_{\alpha 0})$ and $\alpha_{ck}\ \sim\ \mathcal{N}(\alpha_{ck-1},\tau_{c\alpha})$, $\forall k>0$.

### Poll and House Effects

In order to partially pool information among the various scenarios of the same poll, I include a candidate-specific poll effect $(\mu_{cp\[i\]})$, and I also adjust for house effects $(\lambda_{ch\[i\]})$.

### Other Covariates

The vector $X_i$ includes two (standardized) dummies that adjust for the subsample of respondents that the polling firm calculated their estimates on (all respondents, only respondents who are absolutely sure that they will vote in June 2024, or an intermediary subsample). To allow the effect of these covariates to vary as the election date gets closer, these coefficients incorporate a time trend:

$$ \beta_c^{(u)}(date_i)\ =\ \beta_{c0}^{(u)}\ +\ \nu_c(date_i\ -\ 1). $$
