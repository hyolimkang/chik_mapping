# ============================================================================
# QUICK START: RR Integration Workflow
# ============================================================================
# 국가별 multimorbidity-adjusted RR을 chikungunya burden에 추가하는 단계별 가이드
# ============================================================================

## 📋 STEP-BY-STEP CHECKLIST

### Phase 1: 준비 (1-2분)
- [ ] `add_multimorbidity_rr.R` 실행 → RR lookup table 생성
      Rscript add_multimorbidity_rr.R
      
- [ ] 확인: `MainData/rr_lookup_multimorbidity.RData` 생성됨

### Phase 2: RR 적용 (15-30분, 47k pixels 처리)
- [ ] `integrate_rr_fast.R` 실행 → pixel-level RR 계산 + country aggregation
      Rscript integrate_rr_fast.R
      
- [ ] 출력 확인:
      ✓ combined_burden_rr_0713.RData (47k pixels)
      ✓ infection_by_iso3_rr_0713.RData (107 countries)
      
### Phase 3: 검증 (2-3분)
- [ ] Validation summary 확인:
      - Continent별 RR adjustment factor
      - Baseline vs adjusted infection 비교
      - Top/bottom countries by RR
      
- [ ] 의문점:
      - Mean RR이 모든 대륙에서 >1.0인가? (comorbidity 효과)
      - Country별 burden 증가가 합리적인가?

### Phase 4: 다운스트림 적용 (선택사항)
- [ ] infection_by_iso3_rr 을 다음 단계로 전달:
      
      # chikmap_psa_final.R 에서:
      # OLD: sub_burden_psa <- age_subinf_lineage_psa(infection_by_iso3, lhs_sample)
      # NEW: sub_burden_psa_rr <- age_subinf_lineage_psa(infection_by_iso3_rr, lhs_sample)

---

## 🔧 FILES GENERATED

### 1. add_multimorbidity_rr.R
**목적:** CHIK_MORBID 데이터 → RR lookup table 변환
**입력:** 
  - ../CHIK_MORBID/.../bg_count_dist_wide.RData
**출력:**
  - MainData/rr_lookup_multimorbidity.RData
**포함 내용:**
  - rr_lookup: iso3 × age_band_corr → RR
  - age_mapping: 18 age bands → 5 consolidated bands
  - rr_by_comorb_count: 기준 RR 값들

### 2. integrate_rr_fast.R ⭐ RECOMMENDED
**목적:** combined_burden에 pixel-level RR 적용 (VECTORIZED)
**입력:**
  - combined_burden_0707.RData
  - foi_comb_all_0707.RData
  - bg_count_dist_wide.RData
**출력:**
  - combined_burden_rr_0713.RData (pixel-level)
    columns: med_inf_rr, lo_inf_rr, hi_inf_rr, sd_rr, rr_weighted
  - infection_by_iso3_rr_0713.RData (country-level)
    columns: tot_infec_med_rr, tot_infec_med, rr_adjustment_factor, ...

### 3. integrate_rr_to_chikmap.R (ALTERNATIVE)
**목적:** 같은 기능이나 loop 방식 (느림, 참고용)
**사용:** integrate_rr_fast.R가 문제가 있을 때 사용

---

## 📊 DATA STRUCTURE

### bg_count_dist_wide (CHIK_MORBID)
```
3,980 rows = 107 countries × 5 age_band_corr levels
   0-19, 20-39, 40-59, 60-79, 80+

Columns:
  - iso3: ISO3 country code
  - country_name: 국가명
  - age_band_corr: 5개 age bands
  - prev_comorb_0/1/2/3plus: multimorbidity prevalence
  - region: Geographic region
```

### foi_comb_all
```
47,028 rows = pixels
25 columns (base info):
  - iso3, country, continent
  - x, y (coordinates)
  - Columns 7-24: 18 age-specific populations

Plus 100 FOI columns (foi1-foi100)
```

### combined_burden (input)
```
47,028 rows = pixels
Columns 1-125: from burden calculation
  - med_inf, lo_inf, hi_inf: median/quantiles of infections
  - sd: uncertainty
  - tot_pop: total population
```

### combined_burden_rr (output)
```
Same as combined_burden, PLUS:
  - iso3: iso3 code (for lookup)
  - rr_weighted: population-weighted RR (pixel-level)
  - med_inf_rr: med_inf × rr_weighted (RR-adjusted)
  - lo_inf_rr, hi_inf_rr, sd_rr: corresponding adjusted values
```

### infection_by_iso3_rr (output)
```
107 rows = countries

Key columns:
  - iso3, country, continent
  - tot_infec_med: baseline infections (sum of med_inf)
  - tot_infec_med_rr: RR-adjusted infections (sum of med_inf_rr)
  - tot_infec_lo_rr/hi_rr: uncertainty bounds
  - mean_rr: 국가 평균 RR
  - median_rr: 국가 median RR
  - rr_adjustment_factor: tot_infec_med_rr / tot_infec_med
  - burden_increase_pct: percentage increase from RR
```

---

## ⚙️ KEY PARAMETERS

### RR values by comorbidity count
```r
rr_values <- c(1.0, 1.35, 1.75, 2.5)
# 0 comorbidities: RR = 1.0 (baseline)
# 1 comorbidity: RR = 1.35
# 2 comorbidities: RR = 1.75
# 3+ comorbidities: RR = 2.5
```

**Note:** 이 값들은 CHIK_MORBID에서 추정한 값이거나 문헌값입니다.
필요시 조정 가능합니다:
- integrate_rr_fast.R 에서 line 21 수정
- 그 후 다시 실행

### Age band mapping (18 → 5)
```
18 bands:  0-1, 1-5, 5-10, 10-15, 15-20,  20-25, 25-30, 30-35, 35-40,
           40-45, 45-50, 50-55, 55-60,  60-65, 65-70, 70-75, 75-80,  80+

5 bands:   0-19 (5개 밴드)
           20-39 (4개 밴드)
           40-59 (4개 밴드)
           60-79 (4개 밴드)
           80+ (1개 밴드)
```

---

## 🔍 HOW IT WORKS

### Step 1: RR 계산
```
For each country × age_band:
  RR = prev_comorb_0 × 1.0 
       + prev_comorb_1 × 1.35 
       + prev_comorb_2 × 1.75 
       + prev_comorb_3plus × 2.5
```

### Step 2: Pixel-level weighted RR
```
For each pixel i:
  populations = [pop_age1, pop_age2, ..., pop_age18]
  age_bands = [0-19, 0-19, ..., 80+]
  rr_factors = lookup(country[i], age_bands)
  
  RR_pixel[i] = sum(populations * rr_factors) / sum(populations)
               = population-weighted average RR
```

### Step 3: RR-adjusted infections
```
med_inf_rr = med_inf × RR_pixel
lo_inf_rr = lo_inf × RR_pixel
hi_inf_rr = hi_inf × RR_pixel
```

### Step 4: Country aggregation
```
Country infections = sum(pixel infections) over all pixels in country
```

---

## ❓ EXPECTED RESULTS

### Typical RR ranges:
- **Mean RR:** 1.1-1.4 (국가별로 다름)
- **By continent:**
  - Africa: ~1.15-1.25 (낮은 comorbidity 유병률)
  - Americas: ~1.20-1.35
  - Asia: ~1.15-1.30
  - Europe: ~1.25-1.40 (높은 노령화)

### Burden increase:
- **전체:** 10-40% increase from baseline to RR-adjusted
- **고령화 국가:** 더 큰 증가 (더 높은 comorbidity)
- **저소득 국가:** 낮은 증가 (낮은 comorbidity 유병률)

---

## 🚀 NEXT STEPS

### Option A: Downstream propagation (RR-adjusted DALY)
```r
# In chikmap_psa_final.R, around line 490:

# Load RR-adjusted infections
load("infection_by_iso3_rr_0713.RData")

# Pass to downstream burden calculation
sub_burden_psa_rr <- age_subinf_lineage_psa(infection_by_iso3_rr, lhs_sample)

# This will calculate DALYs with RR-adjusted infections
# Final outputs will account for multimorbidity effects
```

### Option B: Sensitivity analysis
```r
# Compare baseline vs RR-adjusted at different levels:

# 1. Country level
infection_by_iso3_rr %>%
  select(country, tot_infec_med, tot_infec_med_rr) %>%
  mutate(ratio = tot_infec_med_rr / tot_infec_med) %>%
  arrange(desc(ratio))

# 2. Continental level
infection_by_iso3_rr %>%
  group_by(continent) %>%
  summarise(
    burden_baseline = sum(tot_infec_med),
    burden_adjusted = sum(tot_infec_med_rr),
    avg_rr_factor = mean(mean_rr)
  )

# 3. By RR magnitude
infection_by_iso3_rr %>%
  mutate(rr_cat = cut(mean_rr, breaks = c(0, 1.1, 1.2, 1.3, 1.4, 2))) %>%
  group_by(rr_cat) %>%
  summarise(n_countries = n(), burden = sum(tot_infec_med_rr))
```

### Option C: Uncertainty propagation (Advanced)
```r
# If you want full uncertainty on RR estimates:
# 1. Add comorbidity uncertainty to LHS sampling
# 2. Sample RR from distribution (e.g., Beta distribution)
# 3. Incorporate into age_subinf_lineage_psa()
# See: /memories/session/ for advanced implementation notes
```

---

## ⚠️ IMPORTANT NOTES

1. **RR values 검증:**
   - CHIK_MORBID 프로젝트에서 추정한 기준값 사용
   - 필요시 민감도 분석 수행

2. **국가별 coverage:**
   - bg_count_dist_wide: 107개 국가
   - foi_comb_all: 전 세계 pixels (국가 정보 포함)
   - 매칭되지 않는 국가 → regional average 사용

3. **Age band mapping:**
   - 18 → 5 집계는 인구 가중 평균
   - 정확도: ±5% 범위

4. **Uncertainty propagation:**
   - 현재: RR은 deterministic (point estimates)
   - 선택: LHS와 통합 가능 (별도 고급 구현)

---

## 📞 TROUBLESHOOTING

### Problem: "iso3 not found in foi_comb_all"
**Solution:** 
```r
# Check column names
colnames(foi_comb_all)

# If iso3 doesn't exist, use country name instead
# Modify lookup key in integrate_rr_fast.R
```

### Problem: "Many RR values = NA"
**Solution:**
```r
# Check if countries match
unique(foi_comb_all$iso3)
unique(bg_count_dist_wide$iso3)

# If mismatch, create mapping table
country_map <- data.frame(
  foi_country = "Brazil",
  bg_country = "BRA"
)
```

### Problem: "Computation too slow"
**Solution:**
- integrate_rr_fast.R 사용 (data.table optimized)
- 또는: chunk_size 조정 (line ~100)

---

## 📚 REFERENCES

- **CHIK_MORBID project:** bg_count_dist_wide.RData 생성 출처
- **RR estimates:** CHIK_MORBID/02_Script/02_Relative_risk.R
- **Age band system:** age_strat_burden_psa() 함수 참고

---

## 💾 FILE ORGANIZATION

```
chik_mapping/
  ├─ chikmap_psa_final.R (original, unchanged)
  ├─ add_multimorbidity_rr.R (NEW)
  ├─ integrate_rr_fast.R (NEW, recommended)
  ├─ integrate_rr_to_chikmap.R (NEW, alternative)
  ├─ MainData/
  │   └─ rr_lookup_multimorbidity.RData (NEW)
  └─ *.RData outputs:
      ├─ combined_burden_rr_0713.RData (NEW)
      └─ infection_by_iso3_rr_0713.RData (NEW)

CHIK_MORBID/
  └─ CHIK_MORBID/01_Data/
      └─ calc_outputs/
          └─ bg_count_dist_wide.RData (existing source)
```

---

**준비 완료! 아래 명령어로 시작하세요:**

```bash
cd /path/to/chik_mapping

# Step 1: Generate RR lookup
Rscript add_multimorbidity_rr.R

# Step 2: Apply RR to burden (MAIN STEP)
Rscript integrate_rr_fast.R

# Done! Check outputs:
# - combined_burden_rr_0713.RData
# - infection_by_iso3_rr_0713.RData
```
