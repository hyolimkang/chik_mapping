# Assumptions

This document has two parts.

Part I (Sections 1–14) sets out the overall analytical framework for a
country-level benefit–risk analysis of chikungunya vaccination, comparing
expected health benefits against vaccine-attributable risks across
countries, age groups, and vaccination strategies (routine, traveller, and
outbreak-response immunisation).

Part II (Section 8, subsections 8.1–8.14) documents the comorbidity-adjusted
hospitalisation burden model that is currently implemented in this
repository. That model supplies the age- and comorbidity-specific
hospitalisation risk used inside the "common clinical outcome model" step
of the Part I framework (Section 8), so its assumptions are nested there
rather than kept as a separate document.

**Note on notation.** Sections 1–7 and 9–14 index country by $c$, age group
by $a$, and comorbidity/multimorbidity category by $m$. Sections 8.1–8.14
retain the notation of the original comorbidity burden analysis, in which
country is indexed by $j$ and comorbidity count is indexed by $c$. These
are equivalent (country: $c \equiv j$; comorbidity category: $m \equiv c$
within 8.1–8.14) — they are not merged here to avoid introducing transcription
errors into already-validated formulas.

---

## 1. Overall analytical framework

The objective is to compare the expected health benefits of chikungunya
vaccination with vaccine-attributable risks across countries, age groups,
and vaccination strategies. The principal benefits will be expressed as
symptomatic cases, hospitalisations, deaths, chronic outcomes, and
disability-adjusted life years (DALYs) averted. Vaccine risks will include
vaccine-attributable serious adverse events, deaths, and corresponding
DALYs.

For country $c$, age group $a$, and vaccination strategy $s$, the general
structure will be:

$$
\text{Benefit}_{c,a,s} = V_{c,a,s} \times P(\text{infection during the relevant risk period})_{c,a,s} \times VE_s \times P(\text{outcome}\mid\text{infection})_a
$$

where $V_{c,a,s}$ is the number vaccinated and $VE_s$ represents the
relevant vaccine effect against infection, symptomatic disease, or severe
disease.

Vaccine-attributable risk will be estimated as:

$$
\text{Vaccine risk}_{c,a,s} = V_{c,a,s} \times P(\text{vaccine-attributable adverse outcome})_a
$$

Benefit–risk ratios will be calculated separately for hospitalisations,
deaths, and DALYs:

$$
BRR_{c,a,s} = \frac{\text{health outcomes averted}_{c,a,s}}{\text{vaccine-attributable outcomes}_{c,a,s}}
$$

Uncertainty in epidemiological, clinical, vaccine-effectiveness, and
vaccine-safety parameters will be propagated through probabilistic
sensitivity analysis. Results will include median benefit–risk ratios,
uncertainty intervals, and the probability that the benefit–risk ratio
exceeds one.

---

## 2. Role of the catalytic model

The existing catalytic model will be retained to estimate:

1. Country-level force of infection;
2. Age-specific baseline susceptibility or prior immunity;
3. The proportion of vaccine recipients who remain susceptible at the time
   of vaccination; and
4. Long-term cumulative infection risk for routine vaccination.

However, the catalytic model will not be used as the sole generator of
age-specific infection burden for all vaccination strategies.

Under a constant force of infection, the instantaneous infection rate at
age $a$ is proportional to:

$$
\lambda_c \exp(-\lambda_c a)
$$

This necessarily produces declining infection rates with increasing age
because older individuals have had more time to acquire immunity. Although
this pattern is internally consistent with a simple endemic catalytic
model, it may not adequately represent chikungunya transmission in
settings characterised by episodic introductions, long inter-epidemic
periods, or explosive outbreaks.

The catalytic model will therefore primarily provide the baseline
susceptible fraction:

$$
S_{c,a} = \exp(-\lambda_c a)
$$

Strategy-specific models will then estimate infection risk during the
relevant vaccination and exposure period.

---

## 3. Strategy-specific infection-risk models

### 3.1 Routine vaccination

Routine vaccination will be evaluated using a cohort-based model rather
than by applying the current annual age-specific infection estimates
directly.

For vaccination at age $a_v$, the proportion susceptible at vaccination
will be estimated from the catalytic model:

$$
S_c(a_v) = \exp(-\lambda_c a_v)
$$

Assuming a constant future force of infection over a time horizon $T$, the
cumulative probability of infection among susceptible vaccine recipients
will be:

$$
P(\text{infection over }T) = 1 - \exp(-\lambda_c T)
$$

The vaccine-preventable infection risk will therefore be:

$$
P(\text{vaccine-preventable infection}) = S_c(a_v) \left[ 1-\exp(-\lambda_c T) \right]
$$

Routine vaccination scenarios will vary:

- Vaccination age;
- Time horizon, such as 5 years, 10 years, or lifetime;
- Vaccine coverage;
- Protection against infection and/or symptomatic disease;
- Duration of protection and waning;
- Baseline serostatus;
- Age-specific disease severity;
- Vaccine-attributable adverse-event risk.

The cohort will be followed over time with ageing, mortality, and, where
relevant, waning vaccine protection.

---

### 3.2 Vaccination of travellers

Traveller vaccination will be evaluated using a duration-specific
infection-risk model. The resident-population age distribution of
infections generated by the catalytic model will not be used directly.

For travel to country $c$ during period $t$, the infection probability
over a trip of duration $d$ will be:

$$
P(\text{travel-associated infection}) = 1-\exp(-\lambda_{c,t}d)
$$

where $\lambda_{c,t}$ represents the infection hazard during the travel
period.

Scenarios will vary:

- Destination-country transmission intensity;
- Travel duration, for example 7, 14, 30, or 90 days;
- Endemic versus outbreak periods;
- Traveller age;
- Baseline immunity;
- Vaccine protection;
- Age-specific risk of severe disease;
- Vaccine-attributable adverse-event risk.

Where reliable time-specific infection hazards are unavailable, destination
risk will be represented using low-, moderate-, and high-transmission
scenarios.

---

### 3.3 Outbreak-response immunisation

Country-specific fitting of weekly surveillance data using a full SEIR
model is not feasible for a global analysis because most countries do not
have sufficiently complete or comparable time-series data.

A simplified semi-mechanistic outbreak-response model will therefore be
used for the primary global analysis. Detailed dynamic transmission models
will be retained for selected data-rich settings and used for calibration
and validation.

---

## 4. Simplified global model for outbreak-response immunisation

### 4.1 Baseline outbreak burden

For country $c$ and age group $a$, infections in the absence of
vaccination will be calculated as:

$$
I^0_{c,a} = N_{c,a} \times S_{c,a} \times AR_c
$$

where:

- $N_{c,a}$ is the population;
- $S_{c,a}$ is the susceptible proportion before the outbreak;
- $AR_c$ is the attack rate among susceptible individuals.

The distinction between the attack rate in the total population and the
attack rate among susceptible individuals will be maintained explicitly to
avoid applying susceptibility twice.

Because reliable country-specific outbreak attack rates are not available
for all countries, attack rate will initially be represented using
scenarios or probability distributions. Illustrative scenarios may include:

$$
AR_c \in \lbrace 0.10,\ 0.25,\ 0.45 \rbrace
$$

representing low-, moderate-, and high-intensity outbreaks among
susceptible individuals.

---

### 4.2 Vaccination timing and the remaining preventable epidemic fraction

The impact of outbreak-response vaccination depends strongly on how much
transmission remains after vaccination has been implemented and vaccine
protection has developed.

Let:

- $\tau$ denote the delay from outbreak onset to vaccination;
- $d_{\text{immune}}$ denote the delay from vaccination to effective
  protection;
- $F_k(t)$ denote the cumulative proportion of outbreak infections
  occurring by time $t$ under outbreak archetype $k$.

The fraction of the epidemic remaining after protection develops will be:

$$
q_k(\tau) = 1-F_k(\tau+d_{\text{immune}})
$$

The directly preventable infections will then be:

$$
I^{\text{averted}}_{c,a,k} = I^0_{c,a} \times v_{c,a} \times VE_{\text{infection}} \times q_k(\tau)
$$

where $v_{c,a}$ is vaccination coverage.

For vaccines that primarily protect against symptomatic disease rather
than infection, the infection process will remain unchanged and prevented
symptomatic outcomes will be estimated as:

$$
C^{\text{averted}}_{c,a,k} = I^0_{c,a} \times v_{c,a} \times VE_{\text{disease}} \times q_k(\tau) \times P(\text{symptomatic}\mid\text{infection})
$$

The primary global analysis will estimate direct effects only. Potential
indirect effects from reduced transmission will be evaluated in
sensitivity analyses and in selected dynamic-model settings.

---

## 5. Outbreak archetypes

Instead of fitting a separate epidemic model to every country, a limited
number of standardised outbreak trajectories will be used.

Potential outbreak archetypes include:

| Archetype | Epidemiological characteristics |
|---|---|
| Fast outbreak | Rapid epidemic growth, early peak, and short duration |
| Intermediate outbreak | Moderate growth rate and duration |
| Slow outbreak | Gradual growth, later peak, and longer duration |
| Multi-wave outbreak | Two or more temporally separated waves |

Each archetype will be characterised by:

- Final attack rate;
- Outbreak duration;
- Time to epidemic peak;
- Time to 25%, 50%, and 75% of cumulative infections;
- The cumulative outbreak function $F_k(t)$;
- The preventable epidemic fraction under different vaccination delays.

Where the appropriate outbreak archetype for a country is uncertain,
outcomes will be averaged across archetypes:

$$
E[I^{\text{averted}}_c] = \sum_k P(k\mid c) I^{\text{averted}}_{c,k}
$$

In the initial analysis, equal weights or broad regional weights may be
used. More informative weights could subsequently be derived from climate,
population density, vector suitability, or historical outbreak
characteristics.

---

## 6. Use of the existing Brazil transmission model

The existing weekly SEIR analysis for Brazil will not be directly
extrapolated to every country. Instead, it will be used as a calibration
dataset for the simplified global model.

Posterior epidemic trajectories from Brazil will be used to estimate:

- Distribution of outbreak attack rates;
- Outbreak duration;
- Time to peak incidence;
- Cumulative infection curves;
- Time to 25%, 50%, and 75% of total infections;
- The proportion of infections remaining under different vaccination
  delays;
- The relationship between vaccination timing and preventable burden;
- The magnitude of indirect effects under infection-blocking vaccination.

For each posterior epidemic trajectory, the preventable fraction will be
calculated for a range of vaccination delays:

$$
q(\tau) = 1-F(\tau+d_{\text{immune}})
$$

This will generate a lookup table or posterior distribution linking the
timing of outbreak-response vaccination to the proportion of the epidemic
that remains preventable.

The resulting simplified timing functions can then be applied to countries
without weekly surveillance data.

---

## 7. Age allocation of infections

The primary outbreak model will not automatically impose the declining
age-specific infection profile generated by the constant-FOI catalytic
model.

Three alternative assumptions will be considered.

### 7.1 Population-proportional exposure

All age groups experience the same per-capita infection risk:

$$
w_{c,a} = \frac{N_{c,a}}{\sum_a N_{c,a}}
$$

### 7.2 Susceptibility-adjusted exposure

Infections are allocated according to both population size and baseline
susceptibility:

$$
w_{c,a} = \frac{N_{c,a}S_{c,a}}{\sum_a N_{c,a}S_{c,a}}
$$

This will be the preferred base case for outbreak-response vaccination.

### 7.3 Catalytic age allocation

The existing catalytic allocation will be retained as a sensitivity
analysis:

$$
w_{c,a} \propto N_{c,a} \lambda_c S_{c,a}
$$

Comparing these assumptions will quantify the extent to which
country-level and age-specific benefit–risk estimates depend on the
catalytic model's declining age profile.

---

## 8. Common clinical outcome model

After estimating infections under each strategy, a common clinical
progression model will be applied:

$$
\text{Infections} \rightarrow \text{Symptomatic cases} \rightarrow \text{Hospitalisations, deaths, and chronic outcomes} \rightarrow \text{DALYs}
$$

For each country, age group, and comorbidity group:

$$
\text{Symptomatic}_{c,a,m} = I_{c,a,m} \times P(\text{symptomatic}\mid\text{infection})
$$

$$
\text{Hospitalised}_{c,a,m} = \text{Symptomatic}_{c,a,m} \times P(\text{hospitalisation}\mid\text{symptomatic},a,m)
$$

where $m$ denotes comorbidity category.

**The comorbidity-adjusted hospitalisation model currently implemented in
this repository provides the $P(\text{hospitalisation}\mid\text{symptomatic},a,m)$
term used at this stage.** Its full derivation, data sources, calibration,
and validation checks are documented in Sections 8.1–8.14 below. Similar
outcome models will be used for mortality and chronic chikungunya where
suitable evidence is available, following the same general structure once
that evidence exists.

### 8.1 Analytical objective

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

### 8.2 Comorbidity-count prevalence

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

### 8.3 Allocation of infections by comorbidity count

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

### 8.4 Symptomatic infections

Symptomatic cases are calculated as:

$$
Y_{j,a,c}=I_{j,a,c}s_r,
$$

where $s_r$ is the symptomatic fraction for region $r$.

In the current implementation, the symptomatic fraction varies by region
or continent but not by age or comorbidity count.

The implemented assumption is therefore:

$$
P(Y=1\mid I=1,r,a,c) = P(Y=1\mid I=1,r).
$$

In particular, comorbidity is assumed not to affect the probability of
developing symptomatic disease after infection.

If age-specific symptomatic fractions are intended, they must be
introduced explicitly because they are not currently included in the
burden function.

### 8.5 Estimation of hospitalisation risk ratios

Age-specific hospitalisation risk ratios are estimated using
individual-level Brazil SINAN data.

A penalised binomial generalised additive model is fitted to estimate:

$$
P(H=1\mid age,c,sex,S=1),
$$

where $S=1$ indicates inclusion as a reported SINAN case.

The model allows the age-risk relationship to vary by comorbidity count:

$$
\text{logit}\lbrace P(H=1) \rbrace = \alpha_c + f_c(age) + \beta_{\mathrm{sex}}sex.
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
RR_{a,c} = \frac{\bar p^{SINAN}_{a,c}}{\bar p^{SINAN}_{a,0}}.
$$

The no-comorbidity group is the reference and therefore has
$RR_{a,0}=1$ in every age group.

The estimated risk ratios vary by both age and comorbidity count. They are
not country-specific.

### 8.6 Transportability of the SINAN risk ratios

The SINAN analysis estimates:

$$
RR^{SINAN}_{a,c} = \frac{P(H=1\mid S=1,a,c)}{P(H=1\mid S=1,a,c=0)}.
$$

The burden model requires:

$$
RR^{target}_{a,c} = \frac{P(H=1\mid symptomatic,a,c)}{P(H=1\mid symptomatic,a,c=0)}.
$$

The primary transportability assumption is:

$$
RR^{SINAN}_{a,c} \approx RR^{target}_{a,c}.
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

### 8.7 External marginal hospitalisation probability

Let:

$$
m_a = P(H=1\mid symptomatic,a)
$$

denote the externally derived age-specific marginal hospitalisation
probability.

This probability is marginal over comorbidity status and is treated as
the absolute age-specific hospitalisation risk that the model must
preserve.

The external marginal probability determines the overall absolute burden,
whereas the Brazil SINAN risk ratios determine how that burden is
distributed across comorbidity groups.

### 8.8 Weighted-RR calibration

For country $j$ and age group $a$, the prevalence-weighted risk ratio is:

$$
W_{j,a} = \sum_c \pi_{j,a,c}RR_{a,c}.
$$

The implied hospitalisation probability for the no-comorbidity reference
group is:

$$
p_{j,a,0} = \frac{m_a}{W_{j,a}}.
$$

The hospitalisation probability for each comorbidity group is:

$$
p_{j,a,c} = p_{j,a,0}RR_{a,c}.
$$

This calibration ensures that:

$$
\sum_c \pi_{j,a,c}p_{j,a,c} = m_a.
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

### 8.9 Hospitalisation burden

Hospitalised cases within each comorbidity group are calculated as:

$$
H_{j,a,c} = Y_{j,a,c}p_{j,a,c}.
$$

The total number of hospitalised cases in an age group is:

$$
H_{j,a} = \sum_c H_{j,a,c}.
$$

Because the weighted-RR calibration preserves the marginal hospitalisation
probability, this total should agree with the burden obtained without
stratifying by comorbidity, subject to numerical precision.

### 8.10 Composition of hospitalisation burden

For the stacked burden figure, each comorbidity group's contribution is
calculated using the same total age-group population as the denominator:

$$
B_{j,a,c} = \frac{H_{j,a,c}}{N_{j,a}} \times10{,}000.
$$

The total height of the stacked bar is:

$$
B_{j,a} = \sum_c B_{j,a,c} = \frac{H_{j,a}}{N_{j,a}} \times10{,}000.
$$

Thus, the total bar height represents the overall hospitalisation burden
per 10,000 people in that age group, and each coloured segment represents
the burden occurring within one comorbidity-count group.

These segments describe burden composition, not causal attribution.

### 8.11 Causal excess burden

The current stacked figure should not be interpreted as the number of
hospitalisations caused by comorbidity.

A causal excess-burden calculation would instead compare the observed
comorbidity-specific probability with a counterfactual probability under
no comorbidity:

$$
H^{excess}_{j,a,c} = Y_{j,a,c} \left( p_{j,a,c}-p_{j,a,0} \right).
$$

Such an interpretation would require stronger causal assumptions,
including adequate control of confounding between comorbidity and
hospitalisation severity.

### 8.12 Uncertainty propagation

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
analyses. Once the benefit–risk framework in Sections 1–14 is implemented,
these same unpropagated sources of uncertainty will also limit the
comorbidity-conditional component of Section 8, and should be reported
alongside the framework's other uncertainty sources (Section 1).

### 8.13 Required validation checks

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

### 8.14 Interpretation

The age-specific RR figure describes differences in hospitalisation risk
conditional on age and reported comorbidity count.

The stacked burden figure describes where hospitalisation burden occurs
in the population.

A high risk ratio does not necessarily imply a large population burden.
Population burden depends jointly on the risk ratio, comorbidity
prevalence, infection burden, symptomatic fraction, and the externally
calibrated absolute hospitalisation probability.

---

## 9. Vaccine effects

At least two vaccine-effect mechanisms will be evaluated.

### 9.1 Disease-only protection

The vaccine reduces symptomatic or severe disease but does not prevent
infection or onward transmission.

Under this mechanism:

- Total infections remain unchanged;
- Symptomatic cases, hospitalisations, deaths, and DALYs are reduced among
  vaccinated individuals;
- No indirect protection is included.

### 9.2 Infection-blocking protection

The vaccine prevents infection and may consequently reduce onward
transmission.

The primary global model will include direct infection prevention only.
Indirect effects will be:

1. Excluded from the conservative base case;
2. Included through scenario-based multipliers in sensitivity analysis; and
3. Estimated dynamically in the Brazil SEIR model and other selected
   data-rich settings.

This distinction is important because the simplified static model can
estimate direct effects reliably, whereas indirect effects depend on
epidemic timing, effective reproduction number, coverage, and spatial
transmission.

---

## 10. Vaccine-attributable risk

Vaccine-attributable adverse outcomes will be estimated by age group:

$$
AE_{c,a,s} = V_{c,a,s} \times P(AE\mid a)
$$

Outcomes will include:

- Serious adverse events;
- Vaccine-attributable deaths;
- DALYs attributable to vaccine adverse events.

Where safety evidence differs substantially by age, separate risk
estimates will be used for younger and older adults.

Benefit and vaccine risk will always be calculated using the same number
of vaccine recipients.

---

## 11. Main scenarios

The global outbreak-response analysis will vary the following parameters:

| Component | Example scenarios |
|---|---|
| Attack rate among susceptible individuals | Low, moderate, high |
| Outbreak archetype | Fast, intermediate, slow, multi-wave |
| Vaccination start | Before onset, onset +2, +4, +6, or +8 weeks |
| Time to protection | Immediate, 7 days, 14 days, or evidence-based value |
| Coverage | 10%, 30%, 50%, 70%, or 90% |
| Rollout duration | Immediate, 2 weeks, 4 weeks, or 8 weeks |
| Vaccine mechanism | Disease-only or infection-blocking |
| Vaccine effectiveness | Posterior or scenario distribution |
| Age eligibility | Adults, older adults excluded, or alternative age ranges |
| Age allocation | Population-proportional, susceptibility-adjusted, catalytic |
| Indirect effects | None, modest, or dynamically estimated |
| Initial immunity | Posterior distribution from the catalytic model |

---

## 12. Primary outputs

Results will be generated by country, age group, comorbidity group, and
vaccination strategy.

Primary outputs will include:

1. Infections averted;
2. Symptomatic cases averted;
3. Hospitalisations averted;
4. Deaths averted;
5. Chronic outcomes averted;
6. DALYs averted;
7. Vaccine-attributable serious adverse events;
8. Vaccine-attributable deaths;
9. Vaccine-attributable DALYs;
10. Benefit–risk ratios;
11. DALYs averted per 10,000 vaccinated;
12. Probability that the benefit–risk ratio exceeds one;
13. Minimum outbreak attack rate required for a favourable benefit–risk
    profile;
14. Maximum vaccination delay consistent with a favourable benefit–risk
    profile.

For outbreak-response immunisation, a particularly useful decision metric
will be the maximum feasible response delay:

$$
\tau^*_{c,a} = \max \left\lbrace \tau : BRR_{c,a,\tau}>1 \right\rbrace
$$

This represents the latest time at which vaccination can begin while
retaining a favourable benefit–risk profile under a specified outbreak
scenario.

---

## 13. Tiered modelling strategy

The analysis will use a tiered modelling framework.

### Tier 1: Global screening analysis

A simplified model will be applied consistently to all countries:

$$
\text{Population} \rightarrow \text{Susceptibility} \rightarrow \text{Outbreak attack rate} \rightarrow \text{Vaccination timing} \rightarrow \text{Outcomes averted} \rightarrow \text{Benefit–risk}
$$

This tier is intended for cross-country comparison and policy screening.

### Tier 2: Dynamic calibration and validation

Detailed SEIR or renewal models will be used in selected data-rich
settings, initially including Brazil, to estimate:

- Epidemic trajectories;
- Vaccination timing effects;
- Indirect protection;
- Rollout effects;
- Accuracy of the simplified global model.

Predictions from Tier 1 will be compared with the corresponding
dynamic-model estimates.

### Tier 3: Country-specific refinement

Where additional country-level surveillance, serological, or outbreak data
are available, the global screening model may be updated using
country-specific attack rates, epidemic curves, reporting rates, or
transmission parameters.

---

## 14. Initial implementation plan

The analysis will be implemented sequentially.

### Step 1: Construct a common country–age input dataset

The master dataset will contain:

- Population;
- Force of infection;
- Baseline susceptibility;
- Comorbidity prevalence;
- Symptomatic fraction;
- Hospitalisation risk;
- Mortality risk;
- Chronic-disease risk;
- DALY weights and durations;
- Vaccine-effectiveness parameters;
- Vaccine-safety parameters.

### Step 2: Separate epidemiological risk modules by strategy

Three functions will be developed:

```r
estimate_routine_risk()
estimate_travel_risk()
estimate_ori_risk()
```
