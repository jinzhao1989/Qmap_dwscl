## Author: jz1989@nwafu.edu.cn
## Date: 2020-03-07

rm(list = ls())
vnms <- "pre"  # ��Ҫ���߶ȵ��������vrnm,�Խ�ˮpreΪ��

# 1.���߶ȿռ䡢ʱ�估������Ϣ
lats <- c(-90,90)   # ���߶�����γ�ȷ�Χ(S,N)
lons <- c(-180,180) # ���߶����򾭶ȷ�Χ(W,E)
yrs_mod <- c(1979,2018) # �۲�(ERA)��ģ��(CRU)���ݷ�Χ1979-2018
yrs_dwn <- c(1901,1978) # ���߶�����(CRU)ʱ�䷶Χ1901-2018
days <- c(1,365)  # ���߶����ڷ�Χ

CRU_pt <- "D:/CRU"  # CRU���ݴ洢·��
ERA_pt <- "D:/ERA"  # ERA���ݴ洢·��

# 2.���߶Ȳ���
# ����R�� (�˴����г��ؼ�R��)
#   qmap      ��λ��ӳ�����
#   splines   ��������
imp <- lapply(c("qmap","splines","raster"), require, character.only=TRUE)
# �����ø߿ռ�ֱ���ERA(0.125��)�������ݴ�������������ѭ������ֹλ��
#   CRU�ֱ���Ϊ0.5��, ERA����Ӧλ�����к���֮���ѭ���м���
rws_era <- c(1,1440) # ERA������ֹ�к�
cls_era <- c(1,2880) # ERA������ʼ�к�
# Qmap��λ�����㲽��
qstp <- 0.01
# ���е��ú���
ncrs <- 20

# 3.�������
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # WGS84��������ϵ����
fsp <- .Platform$file.sep
out_pt <- "D:/Downscale"  # ���·��
bv <- -9999  # ����ֵ

# 4.���л����߶�
# ��ѭ��
for (day in days[1]:days[2]) {
  # ��ȡ��day��CRU���й۲���yrs_mod���ճ߶���������
  CRU_mod_dts <- lapply(yrs_mod[1]:yrs_mod[2], FUN = function(iyr){
    # CRU pre��iyr���day�������ļ��� CRU_pre_iyr_day.tif
    CRU_mod_dt <- raster(CRU_pre_iyr_day.tif)  # ����ݽ�һ����ȡΪmatrix��array
  })  # �˴��ٶ��ѽ�list���ͽ��ת��Ϊarray������ά����
  
  # ��ȡ��day��ERA���й۲���yrs_mod���ճ߶���������
  ERA_obv_dts <- lapply(yrs_mod[1]:yrs_mod[2], FUN = function(iyr){
    # ERA pre��iyr���day�������ļ��� ERA_pre_iyr_day.tif
    ERA_obv_dt <- raster(ERA_pre_iyr_day.tif)  # ͬ��
  })  # ͬ��
  
  # ��ȡ��day��CRU���н��߶���yrs_dwn���ճ߶���������
  CRU_dwn_dts <- lapply(yrs_dwn[1]:yrs_dwn[2], FUN = function(iyr){
    # CRU pre��iyr���day�������ļ��� CRU_pre_iyr_day.tif
    CRU_dwn_dt <- raster(CRU_pre_iyr_day.tif)  # ͬ��
  })  # ͬ��
  
  # �½�����ά����dwns_day�洢��day��yrs_dwn�����꽵�߶Ƚ��
  dwns_day <- array(data = NA,dim = c(diff(rws_era)+1,diff(cls_era)+1),diff(yrs_dwn)+1)
  # ��ѭ��,��õ�day��yrs_dwn�����꽵�߶���ά�������dwn_day
  for (rw_era in rws_era) {
    # ����CRU(0.5��)��Ӧ�к�
    rw_cru <- ceiling(rw_era/(0.5/0.125))
    # ��ѭ��
    for (cl_era in cls_era) {
      # ����CRU(0.5��)��Ӧ�к�
      cl_cru <- ceiling(cl_era/(0.5/0.125))
      # ��ȡָ�����к�λ�ø���������������
      CRU_mod_dt <- CRU_mod_dts[rw_cru,cl_cru,]
      ERA_obv_dt <- ERA_obv_dts[rw_era,cl_era,]
      CRU_dwn_dt <- CRU_dwn_dts[rw_cru,cl_cru,]
      # ����Qmap������ù۲�ʱ���CRU��ERA��Ӧ��ϵ
      try_q <- try(fit_qm <- fitQmapDIST(ERA_obv_dt,CRU_mod_dt, qstep=qstp),silent = T)
      if (class(try_q)=="try-error") {
        # ���Qmap����������,���滻Ϊ����ƽ��������Ϲ�ϵ(cubic smoothing spline)
        try_sp <- try(fit_sp <- smooth.spline(CRU_mod_dt,ERA_obv_dt),silent = T)
        if (class(try_sp)=="try-error") {
          # �������2�ַ�����������,�����NAֵ
          down_data <- rep(NA,diff(yrs_dwn)+1)
        } else {
          # �������ƽ��������������,���ø÷�����Ͻ��߶�
          down_data <- predict(fit_sp,CRU_dwn_dt)$y
        }
      } else {
        # ���Qmap��������,����Qmap��õĹ�ϵ���߶�
        down_data <- doQmapDIST(CRU_dwn_dt, fit_qm)
      }
      dwns_day[rw_era,cl_era,] <- down_data
    }
  }
  # ���dwns_dayÿ��������
  for (iyr in yrs_dwn[1]:yrs_dwn[2]) {
    dwn_day1 <- dwns_day[,,(iyr-yrs_dwn[1]+1)]
    dwn_day1[is.na(dwn_day1)] <- bv
    dwn_day1 <- raster(x = dwn_day1,xmn=lons[1], xmx=lons[2],
                       ymn=lats[1], ymx=lats[2],crs=wgs84)
    writeRaster(x = dwn_day1,filename = 
                  paste0(out_pt,fsp,vnms,"_dwscl_",iyr,sprintf("%03d",day),".tif"))
  }
}

# End