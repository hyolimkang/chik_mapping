# Inspect chikv_fatal_hosp_rate.RData structure
load('MainData/chikv_fatal_hosp_rate.RData')

cat('=== 로드된 객체 ===\n')
print(ls())

cat('\n=== hosp 구조 및 첫 10개 값 ===\n')
print(class(hosp))
print(str(hosp))
print(head(hosp, 10))

cat('\n=== fatal 구조 및 첫 10개 값 ===\n')
print(class(fatal))
print(str(fatal))
print(head(fatal, 10))

cat('\n=== nh_fatal 구조 및 첫 10개 값 ===\n')
print(class(nh_fatal))
print(str(nh_fatal))
print(head(nh_fatal, 10))

cat('\n=== 값의 범위 ===\n')
cat('hosp 범위:', min(hosp), '~', max(hosp), '\n')
cat('fatal 범위:', min(fatal), '~', max(fatal), '\n')
cat('nh_fatal 범위:', min(nh_fatal), '~', max(nh_fatal), '\n')

# Check if there's any metadata
cat('\n=== 환경 내 모든 변수 ===\n')
print(ls(all.names = TRUE))
