# age_strat_burden_psa 함수 분석 및 Multimorbidity 통합 가이드

## 1. age_strat_burden_psa 함수 상세 분석

### 1.1 함수 위치 및 목적
- **파일**: `Functions/BurdenFunctions_v2.R` (line 224~265, NA 처리 버전)
- **목적**: 100개의 FOI 샘플을 25개씩 4개 배치로 나누어 age-stratified infection 계산

### 1.2 입력 데이터
```r
burden_1 <- age_strat_burden_psa(df = foi_comb_all, 26, 50)
```
- `df`: foi_comb_all (픽셀별 FOI + age-specific population)
- `start_col`: 26, `end_col`: 50 (FOI 컬럼 범위)
  - foi1~foi25: burden_1
  - foi26~foi50: burden_2
  - foi51~foi75: burden_3
  - foi76~foi100: burden_4

### 1.3 핵심 계산 프로세스

#### Step 1: Age bands 정의
```r
ages <- c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)
l_lim <- ages
u_lim <- c(1, ages[3:18]-1, 100)
# 결과: 18개 age groups
```

#### Step 2: 각 FOI 샘플별 처리 (j = start_col to end_col)
```r
for(i in 1:nrow(df)) {  # 각 픽셀
  # Age band별 incidence rate 계산
  incidence_rates <- mapply(
    calc_incidence, 
    FOI = rep(df[i, foi_col], length(ages)), 
    l_lim, u_lim
  )
  
  # Age-specific population 추출
  age_group_pop <- as.numeric(df[i, 10:27])  # 18개 age 그룹
  
  # 감염자 수 계산
  infection <- mapply(
    incidence_to_numbs, 
    incidence = incidence_rates, 
    n_j = age_group_pop
  )
  
  # 저장
  infections_per_age_band[i, ] <- infection
  total_infection[i] <- sum(infection)
}
```

#### Step 3: incidence 계산 함수
```r
calc_incidence <- function(FOI, l_lim, u_lim) {
  # Step 1: Age band 내 감염 확률
  p_I <- exp(-FOI * l_lim) - exp(-FOI * u_lim)
  
  # Step 2: Incidence rate (연령당 평균 위험)
  incidence <- p_I / (u_lim - l_lim)
  
  return(incidence)
}
```

### 1.4 출력 데이터
각 요소는 list with 3 components:
```
[[j]]$updated_df           # 원본 df + total_infection 컬럼
[[j]]$infection_per_band   # (nrow(df) × 18) matrix - age별 감염자 수
[[j]]$incidence_per_band   # (nrow(df) × 18) matrix - age별 incidence rate
```

최종 통합:
```r
tot_inf <- do.call(cbind, lapply(all_results, function(df) df$updated_df$total_infection))
# 결과: nrow(df) × 100 matrix (100개 FOI 샘플별 total infections)

combined_burden$med_inf <- apply(tot_inf_cols, 1, median, na.rm = T)
combined_burden$lo_inf <- apply(tot_inf_cols, 1, function(x) quantile(x, probs = 0.025))
combined_burden$hi_inf <- apply(tot_inf_cols, 1, function(x) quantile(x, probs = 0.975))
```

---

## 2. CHIK_MORBID 프로젝트의 Relative Risk 데이터

### 2.1 데이터 구조
**경로**: `CHIK_MORBID/CHIK_MORBID/02_Script/02_Relative_risk.R`

### 2.2 주요 변수
- **Comorbidities** (7가지):
  - Diabetes
  - Hypertension
  - Hepatopathy (간질환)
  - Renal disease
  - Hematologic disease
  - Peptic ulcer
  - Autoimmune disease

- **Multimorbidity groups**:
  - 0: 질환 없음
  - 1: 1개 질환
  - 2: 2개 질환
  - 3+: 3개 이상

- **Outcomes** (severe outcome):
  - Hospitalization = "yes" OR Death = TRUE

### 2.3 상대 위험도(RR) 모델
```r
# 기본 모델: age × comorbidity 상호작용
fit_condition_only <- glm(
  severe ~ ns(age_years, df = 4) + 
    diabetes + hypertension + hepatopathy + renal_disease +
    hematologic + peptic_ulcer + autoimmune + sex,
  data = analysis_df,
  family = binomial()
)

# Multimorbidity count 모델
fit_morbid_burden <- glm(
  severe ~ ns(age_years, df = 4) +
    comorb_count_group +
    sex,
  data = analysis_df,
  family = binomial()
)
```

### 2.4 RR 계산 방식
```r
rr_all_severe <- pred_all_severe |>
  pivot_wider(names_from = condition_status, values_from = p_severe) |>
  mutate(
    rr = yes/no,  # Relative Risk = P(severe|condition) / P(severe|no condition)
    risk_difference = yes - no
  )
```

결과: age × condition 매트릭스 (18 age groups × 7 conditions)

---

## 3. Downstream Burden 계산 전략

### 3.1 현재 구조 (기존)
```
FOI (100 samples) 
    ↓
age_strat_burden_psa() 
    ↓
age_strat_infections (100 × 18 age groups)
    ↓
age_subinf_lineage_psa() 
    ↓
DALY components (YLL, YLD_hosp, YLD_nonhosp, YLD_chronic)
```

### 3.2 제안하는 확장 구조
```
FOI (100 samples) 
    ↓
age_strat_burden_psa() 
    ↓
age_strat_infections (100 × 18 age groups)
    ↓
[NEW] age_strat_infections_with_rr()  ← RR adjustment 추가
      - Multimorbidity count 통합
      - Age-specific RR 적용
      - RR-adjusted infections 계산
    ↓
[MODIFIED] age_subinf_lineage_psa_with_rr()
      - RR-adjusted infections 입력
      - 원래 burden calculation 동일
    ↓
DALY components (RR-adjusted)
```

---

## 4. 구현 방안 (3가지 옵션)

### 옵션 A: Post-hoc RR 적용 (권장)
```r
# Step 1: Population comorbidity prevalence 계산
# CHIK_MORBID의 population prevalence data 활용
# → age × comorbidity_count matrix 생성

# Step 2: RR-adjusted infection rate
# RR_adjusted_infections = base_infections × 
#   (1 + multimorbidity_prevalence × (RR_relative_increase))

# 장점: 기존 코드 최소 수정, modular
# 단점: multimorbidity와 FOI 간 correlation 미고려
```

### 옵션 B: 함수 수정 (Medium complexity)
```r
# burden 계산 함수 수정
age_strat_burden_psa_with_rr <- function(
  df, 
  start_col, 
  end_col,
  rr_matrix,           # age × condition RR matrix
  multimorbidity_df    # pixel-level multimorbidity data
) {
  # 기존 로직 + RR 곱하기
  infection <- infection * rr_weight[age_group]
}
```

### 옵션 C: LHS 통합 (High complexity)
```r
# LHS sample에 multimorbidity 변수 추가
# → 1000개 LHS runs 각각에 age-specific RR 적용
# → downstream burden에서 직접 통합
```

---

## 5. 권장 구현 단계

### Phase 1: 데이터 준비
1. CHIK_MORBID에서 age-specific RR 데이터 추출 및 저장
   ```r
   # Save as: MainData/rr_by_age_condition.RData
   # Structure: data.frame(age, age_band, condition, rr_no, rr_yes, rr_ratio)
   ```

2. 전역 multimorbidity prevalence 추정
   ```r
   # Population-level prevalence by age
   # PNS 2019 Brazil data 활용
   multimorbidity_prev <- read.csv("MainData/multimorbidity_prev_by_age.csv")
   ```

### Phase 2: 함수 개발
```r
# 새 함수 추가 (BurdenFunctions_v2.R)
apply_rr_adjustment <- function(
  infections_age_stratified,  # age_strat_burden_psa 출력
  rr_matrix,
  multimorbidity_prevalence,
  method = "additive"  # or "multiplicative"
) {
  # infection adjustment 로직
}
```

### Phase 3: 통합 및 검증
```r
# 기존 코드에 추가
sub_burden_psa_with_rr <- age_subinf_lineage_psa(
  infection_by_iso3,        # RR-adjusted inputs
  lhs_sample
)
```

---

## 6. 데이터 연결 고리

| Component | Location | Format | Usage |
|-----------|----------|--------|-------|
| Base FOI | chik_mapping/MainData/foi_mat.RData | 100 samples × pixel | age_strat_burden_psa input |
| Age-stratified pop | foi_comb_all[10:27] | 18 age groups per pixel | Denominator |
| RR by age/condition | CHIK_MORBID/02_Script/02_Relative_risk.R | glm predictions | RR adjustment |
| Multimorbidity prev | CHIK_VIM/00_Data/population_comorbidity_pns2019.csv | age × condition | Weight for RR |
| LHS samples | chikmap_psa_final.R line ~150 | 1000 samples | Parameter uncertainty |

---

## 7. 다음 단계

사용자의 의도 및 선호도 확인 필요:
- [ ] Multimorbidity를 deterministic (평균값) vs probabilistic (LHS)로 적용?
- [ ] Country-specific medical care 감소 + RR 증폭 combined?
- [ ] 특정 comorbidity만 focus vs 모두 포함?
- [ ] DALY 컴포넌트 중 어느 것에 RR 적용 (hospitalization, severity, fatality)?
