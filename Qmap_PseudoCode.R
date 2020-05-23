## Author: jz1989@nwafu.edu.cn
## Date: 2020-05-22

rm(list = ls())
vnms <- "pre"  # 需要降尺度的气象变量vrnm,以降水pre为例

# 1.降尺度空间、时间及数据信息
lats <- c(-90,90)   # 降尺度区域纬度范围(S,N)
lons <- c(-180,180) # 降尺度区域经度范围(W,E)
yrs_mod <- c(1979,2018) # 观测(ERA)及模拟(CRU)数据范围1979-2018
yrs_dwn <- c(1901,1978) # 降尺度数据(CRU)时间范围1901-2018
days <- c(1,365)  # 降尺度日期范围

CRU_pt <- "D:/CRU"  # CRU数据存储路径
ERA_pt <- "D:/ERA"  # ERA数据存储路径

# 2.降尺度参数
# 加载R包 (此处仅列出关键R包)
#   qmap      分位数映射后处理
#   splines   样条估计
imp <- lapply(c("qmap","splines","raster"), require, character.only=TRUE)
# 计算获得高空间分辨率ERA(0.125°)格网数据待处理区域行列循环的起止位置
#   CRU分辨率为0.5°, ERA格点对应位置行列号在之后的循环中计算
rws_era <- c(1,1440) # ERA数据起止行号
cls_era <- c(1,2880) # ERA数据起始列号
# Qmap分位数计算步长
qstp <- 0.01

# 3.输出参数
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # WGS84地理坐标系参数
fsp <- .Platform$file.sep
out_pt <- "D:/Downscale"  # 输出路径
bv <- -9999  # 背景值

# 4.并行化降尺度
# 日循环
for (day in days[1]:days[2]) {
  # 读取第day日CRU所有观测年yrs_mod的日尺度网格数据
  CRU_mod_dts <- lapply(yrs_mod[1]:yrs_mod[2], FUN = function(iyr){
    # CRU pre第iyr年第day天数据文件名 CRU_pre_iyr_day.tif
    CRU_mod_dt <- raster(CRU_pre_iyr_day.tif)  # 需根据进一步提取为matrix或array
  })  # 此处假定已将list类型结果转换为array类型三维矩阵
  
  # 读取第day日ERA所有观测年yrs_mod的日尺度网格数据
  ERA_obv_dts <- lapply(yrs_mod[1]:yrs_mod[2], FUN = function(iyr){
    # ERA pre第iyr年第day天数据文件名 ERA_pre_iyr_day.tif
    ERA_obv_dt <- raster(ERA_pre_iyr_day.tif)  # 同上
  })  # 同上
  
  # 读取第day日CRU所有降尺度年yrs_dwn的日尺度网格数据
  CRU_dwn_dts <- lapply(yrs_dwn[1]:yrs_dwn[2], FUN = function(iyr){
    # CRU pre第iyr年第day天数据文件名 CRU_pre_iyr_day.tif
    CRU_dwn_dt <- raster(CRU_pre_iyr_day.tif)  # 同上
  })  # 同上
  
  # 新建空三维数组dwns_day存储第day日yrs_dwn所有年降尺度结果
  dwns_day <- array(data = NA,dim = c(diff(rws_era)+1,diff(cls_era)+1),diff(yrs_dwn)+1)
  # 行循环,获得第day日yrs_dwn所有年降尺度三维结果矩阵dwn_day
  for (rw_era in rws_era) {
    # 计算CRU(0.5°)对应行号
    rw_cru <- ceiling(rw_era/(0.5/0.125))
    # 列循环
    for (cl_era in cls_era) {
      # 计算CRU(0.5°)对应列号
      cl_cru <- ceiling(cl_era/(0.5/0.125))
      # 读取指定行列号位置格点所有年份日数据
      CRU_mod_dt <- CRU_mod_dts[rw_cru,cl_cru,]
      ERA_obv_dt <- ERA_obv_dts[rw_era,cl_era,]
      CRU_dwn_dt <- CRU_dwn_dts[rw_cru,cl_cru,]
      # 尝试Qmap方法获得观测时间段CRU与ERA对应关系
      try_q <- try(fit_qm <- fitQmapDIST(ERA_obv_dt,CRU_mod_dt, qstep=qstp),silent = T)
      if (class(try_q)=="try-error") {
        # 如果Qmap方法不可行,则替换为三次平滑样条拟合关系(cubic smoothing spline)
        try_sp <- try(fit_sp <- smooth.spline(CRU_mod_dt,ERA_obv_dt),silent = T)
        if (class(try_sp)=="try-error") {
          # 如果上述2种方法均不可行,则输出NA值
          down_data <- rep(NA,diff(yrs_dwn)+1)
        } else {
          # 如果三次平滑样条方法可行,则用该方法拟合降尺度
          down_data <- predict(fit_sp,CRU_dwn_dt)$y
        }
      } else {
        # 如果Qmap方法可行,则用Qmap获得的关系降尺度
        down_data <- doQmapDIST(CRU_dwn_dt, fit_qm)
      }
      dwns_day[rw_era,cl_era,] <- down_data
    }
  }
  # 输出dwns_day每层日数据
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
