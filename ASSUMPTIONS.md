# Assumptions

## Analytical objective

The objective is to estimate age-specific chikungunya hospitalisation
burden by comorbidity count and to describe how the total burden is
distributed across people with zero, one, or two or more comorbidities.

The analysis combines:

1. country- and age-specific chikungunya infection estimates;
2. general-population prevalence of comorbidity counts;
3. symptomatic fractions;
4. external age-specific marginal hospitalisation probabilities; and
5. age-specific hospitalisation risk ratios estimated from Brazil SINAN.

The primary analysis describes burden occurring within each comorbidity
group. It does not estimate the causal excess burden attributable to
comorbidity.

## Comorbidity-count prevalence

Country- and age-specific marginal prevalence estimates for individual
conditions are combined with an age-specific correlation matrix using a
Gaussian copula simulation.

For every country and age group, simulated individuals are classified as
having:

- no comorbidity;
- one comorbidity; or
- two or more comorbidities.

The resulting prevalence distribution is denoted by:

$$
\pi_{j,a,c},
$$

where $j$ indexes country, $a$ indexes age group, and $c$ indexes
comorbidity count.

The distributions must satisfy:

$$
\sum_c \pi_{j,a,c}=1.
$$

The copula procedure assumes that the available marginal prevalence
estimates and the selected age-specific correlation matrices adequately
represent the joint distribution of the included conditions.

## Allocation of infections by comorbidity count

Total infections in country $j$ and age group $a$ are allocated to
comorbidity groups according to their prevalence in the general population:

$$
I_{j,a,c}=I_{j,a}\pi_{j,a,c}.
$$

This calculation assumes that, within the same age group, the probability
of chikungunya infection does not differ by comorbidity count:

$$
P(I=1\mid j,a,c)=P(I=1\mid j,a).
$$

Consequently, differences in the number of infections across comorbidity
groups arise from differences in subgroup population size rather than
from differences in infection risk.

This assumption may be violated if comorbidity is associated with exposure,
mobility, housing, mosquito contact, healthcare behaviour, or biological
susceptibility to infection.

## Symptomatic infections

Symptomatic cases are calculated as:

$$
Y_{j,a,c}=I_{j,a,c}s_r,
$$

where $s_r$ is the symptomatic fraction for region $r$.

In the current implementation, the symptomatic fraction varies by region
or continent but not by age or comorbidity count.

The implemented assumption is therefore:

$$
P(Y=1\mid I=1,r,a,c)
=
P(Y=1\mid I=1,r).
$$

In particular, comorbidity is assumed not to affect the probability of
developing symptomatic disease after infection.

If age-specific symptomatic fractions are intended, they must be introduced
explicitly because they are not currently included in the burden function.

## Estimation of hospitalisation risk ratios

Age-specific hospitalisation risk ratios are estimated using individual-level
Brazil SINAN data.

A penalised binomial generalised additive model is fitted to estimate:

$$
P(H=1\mid age,c,sex,S=1),
$$

where $S=1$ indicates inclusion as a reported SINAN case.

The model allows the age-risk relationship to vary by comorbidity count:

$$
\text{logit}\{P(H=1)\}
=
\alpha_c + f_c(age) + \beta_{\mathrm{sex}}sex.
$$

Predicted hospitalisation probabilities are standardised within each
10-year age band by:

1. generating predictions for each single year of age;
2. averaging over the specified sex distribution; and
3. averaging the standardised probabilities across ages within the band.

Let the resulting standardised probability be:

$$
\bar p^{SINAN}_{a,c}.
$$

The age-specific risk ratio is then:

$$
RR_{a,c}
=
\frac{
\bar p^{SINAN}_{a,c}
}{
\bar p^{SINAN}_{a,0}
}.
$$

The no-comorbidity group is the reference and therefore has
$RR_{a,0}=1$ in every age group.

The estimated risk ratios vary by both age and comorbidity count. They are
not country-specific.

## Transportability of the SINAN risk ratios

The SINAN analysis estimates:

$$
RR^{SINAN}_{a,c}
=
\frac{
P(H=1\mid S=1,a,c)
}{
P(H=1\mid S=1,a,c=0)
}.
$$

The burden model requires:

$$
RR^{target}_{a,c}
=
\frac{
P(H=1\mid symptomatic,a,c)
}{
P(H=1\mid symptomatic,a,c=0)
}.
$$

The primary transportability assumption is:

$$
RR^{SINAN}_{a,c}
\approx
RR^{target}_{a,c}.
$$

In words, age-specific hospitalisation risk ratios estimated among
SINAN-reported cases are assumed to apply to all symptomatic chikungunya
infections.

This assumption may be violated if selection into SINAN depends jointly
on comorbidity status and disease severity. For example, people with
comorbidities may seek care or receive diagnostic testing for milder
illness, whereas hospitalised cases may be reported with high probability
regardless of comorbidity.

Such selection could bias the risk ratios in either direction.

Direct correction would require information on ascertainment probabilities
stratified jointly by:

- age;
- comorbidity status; and
- hospitalisation or disease severity.

## External marginal hospitalisation probability

Let:

$$
m_a
=
P(H=1\mid symptomatic,a)
$$

denote the externally derived age-specific marginal hospitalisation
probability.

This probability is marginal over comorbidity status and is treated as
the absolute age-specific hospitalisation risk that the model must
preserve.

The external marginal probability determines the overall absolute burden,
whereas the Brazil SINAN risk ratios determine how that burden is distributed
across comorbidity groups.

## Weighted-RR calibration

For country $j$ and age group $a$, the prevalence-weighted risk ratio is:

$$
W_{j,a}
=
\sum_c
\pi_{j,a,c}RR_{a,c}.
$$

The implied hospitalisation probability for the no-comorbidity reference
group is:

$$
p_{j,a,0}
=
\frac{m_a}{W_{j,a}}.
$$

The hospitalisation probability for each comorbidity group is:

$$
p_{j,a,c}
=
p_{j,a,0}RR_{a,c}.
$$

This calibration ensures that:

$$
\sum_c
\pi_{j,a,c}p_{j,a,c}
=
m_a.
$$

Therefore, separating the population into comorbidity groups does not
change the externally specified marginal hospitalisation probability.

The calibration must be checked to ensure that all estimated subgroup
probabilities remain between zero and one:

$$
0\le p_{j,a,c}\le1.
$$

Values greater than one should be treated as evidence of incompatible
marginal risks, prevalence distributions, or relative risks rather than
being silently truncated.

## Hospitalisation burden

Hospitalised cases within each comorbidity group are calculated as:

$$
H_{j,a,c}
=
Y_{j,a,c}p_{j,a,c}.
$$

The total number of hospitalised cases in an age group is:

$$
H_{j,a}
=
\sum_c H_{j,a,c}.
$$

Because the weighted-RR calibration preserves the marginal hospitalisation
probability, this total should agree with the burden obtained without
stratifying by comorbidity, subject to numerical precision.

## Composition of hospitalisation burden

For the stacked burden figure, each comorbidity group's contribution is
calculated using the same total age-group population as the denominator:

$$
B_{j,a,c}
=
\frac{
H_{j,a,c}
}{
N_{j,a}
}
\times10{,}000.
$$

The total height of the stacked bar is:

$$
B_{j,a}
=
\sum_c B_{j,a,c}
=
\frac{
H_{j,a}
}{
N_{j,a}
}
\times10{,}000.
$$

Thus, the total bar height represents the overall hospitalisation burden
per 10,000 people in that age group, and each coloured segment represents
the burden occurring within one comorbidity-count group.

These segments describe burden composition, not causal attribution.

## Causal excess burden

The current stacked figure should not be interpreted as the number of
hospitalisations caused by comorbidity.

A causal excess-burden calculation would instead compare the observed
comorbidity-specific probability with a counterfactual probability under
no comorbidity:

$$
H^{excess}_{j,a,c}
=
Y_{j,a,c}
\left(
p_{j,a,c}-p_{j,a,0}
\right).
$$

Such an interpretation would require stronger causal assumptions,
including adequate control of confounding between comorbidity and
hospitalisation severity.

## Uncertainty propagation

The uncertainty intervals in the final burden estimates may include
uncertainty in:

- infection estimates;
- symptomatic fractions;
- external marginal hospitalisation probabilities;
- disability weights and disease durations; and
- other Latin hypercube sampling parameters.

If the median SINAN risk ratio is copied unchanged across all simulation
runs, uncertainty in the GAM-derived risk ratios is not propagated into
the final burden intervals.

Similarly, unless explicitly sampled, the following sources of uncertainty
are not propagated:

- uncertainty in individual-condition prevalence estimates;
- uncertainty in the comorbidity correlation matrices;
- Monte Carlo uncertainty from the copula simulation; and
- uncertainty in the transportability of SINAN risk ratios.

These should be reported as limitations or incorporated through sensitivity
analyses.

## Required validation checks

The following numerical checks should be conducted before interpreting the
results:

1. Comorbidity-count prevalence sums to one within every country and age
   group:

$$
\sum_c\pi_{j,a,c}=1.
$$

2. Subgroup infections recover total infections:

$$
\sum_c I_{j,a,c}=I_{j,a}.
$$

3. All hospitalisation probabilities are within zero and one.

4. The prevalence-weighted subgroup probabilities recover the external
   marginal probability:

$$
\sum_c\pi_{j,a,c}p_{j,a,c}=m_a.
$$

5. Subgroup hospitalisations recover the total:

$$
\sum_c H_{j,a,c}=H_{j,a}.
$$

6. The stacked contributions recover the overall hospitalisation rate:

$$
\sum_c B_{j,a,c}=B_{j,a}.
$$

7. All joins between age groups, countries, comorbidity groups, and
   simulation runs are complete and do not produce missing values.

## Interpretation

The age-specific RR figure describes differences in hospitalisation risk
conditional on age and reported comorbidity count.

The stacked burden figure describes where hospitalisation burden occurs
in the population.

A high risk ratio does not necessarily imply a large population burden.
Population burden depends jointly on the risk ratio, comorbidity prevalence,
infection burden, symptomatic fraction, and the externally calibrated
absolute hospitalisation probability.
