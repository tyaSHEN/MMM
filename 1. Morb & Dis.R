library(tidyverse)
library(nnet)
library(doParallel)
library(haven)


# clean data----
# only need to run once, cleaned data will be stored in folder "output"
if(!dir.exists(file.path("output"))){dir.create(file.path("output"))}

# this import is slow so we will convert it to a Rdata file first
# data.Rdata stores a R data.frame named "Raw" converted from "randhrs1992_2018v1_STATA.zip"
# https://hrsdata.isr.umich.edu/data-products/rand-hrs-longitudinal-file-2018

# Raw <- read_dta("randhrs1992_2018v2.dta")
# save(Raw,file="data.RData")

load("data.RData")

# two health dimension ---
ADL = Raw %>% select("rahhidpn",
                     paste0("r",c(2:14),"adla"))

ADL = pivot_longer(ADL,c(2:ncol(ADL)))
ADL$wave= parse_number(ADL$name)
ADL$name = gsub("[0-9]","",ADL$name)
ADL = pivot_wider(ADL,names_from=name,values_from = value)
colnames(ADL)[3]="radl"

DIS = Raw %>% select("rahhidpn",
                     paste0("r",c(2:14),"diabe"),
                     paste0("r",c(2:14),"cancre"),
                     paste0("r",c(2:14),"lunge"),
                     paste0("r",c(2:14),"hearte"),
                     paste0("r",c(2:14),"stroke")) 
DIS = pivot_longer(DIS,c(2:(13*5+1)))
DIS$wave= parse_number(DIS$name)
DIS$name = gsub("[0-9]","",DIS$name)
DIS$value = as.numeric(DIS$value)
DIS = pivot_wider(DIS,names_from=name,values_from = value)
DIS = DIS %>% mutate(rdis = select(., rdiabe:rstroke) %>% rowSums(na.rm = TRUE))
DIS$rdis = ifelse(is.na(DIS$rdiabe) & is.na(DIS$rcancre) & is.na(DIS$rlunge) & is.na(DIS$rhearte) & is.na(DIS$rstroke),NA,DIS$rdis)
DIS = DIS[,c(1,2,8)]


SW = Raw %>% select("rahhidpn",paste0("r",c(2:13),"wtcrnh"))
SW = pivot_longer(SW,c(2:13))
SW$wave= parse_number(SW$name)
SW$name = gsub("[0-9]","",SW$name)
SW = pivot_wider(SW,names_from=name,values_from = value)


dem = Raw %>% select("rahhidpn","raeduc","ragender","raracem","rahispan")
dem$edu = ifelse(dem$raeduc %in% c(2,3),"3",dem$raeduc)
dem$race = dem$raracem
dem$race = ifelse(dem$race == 3, 4 ,dem$race)
dem$race = ifelse(dem$rahispan == 1, 3 ,dem$race)
dem = dem %>% filter(race != 4)

YAW = data.frame(FY = c(1934,1944,1924,1934,1914,1924),
                 FA = c(60,60,70,70,80,80),
                 WA = c(4,9,4,9,4,9))

YAW$LY = YAW$FY+9
YAW$LA = YAW$FA+9                 
save.image(file = "output/dem.Rdata")

for (CO in unique(YAW$FY)) {
  Cohort = Raw %>% select("rahhidpn","rabyear","radyear",paste0("r",c(4:14),"agey_e")) %>% 
    filter(rabyear %in% c(CO:(CO+9))) %>%
    mutate(WD = as.integer(radyear/2-995),WD = ifelse(is.na(WD),14,WD))
  Cohort = pivot_longer(Cohort,c(4:14),values_to = "age")
  Cohort$wave= parse_number(Cohort$name)
  Cohort = Cohort %>% select(-name)
  
  Cohort = left_join(Cohort,ADL)
  Cohort$radl=ifelse(Cohort$wave>Cohort$WD,999,Cohort$radl)
  
  Cohort = left_join(Cohort,DIS)
  Cohort$rdis=ifelse(Cohort$wave>Cohort$WD,999,Cohort$rdis)
  
  Cohort$dis=ifelse(Cohort$rdis <1,"N",ifelse(Cohort$rdis<3,"S",ifelse(Cohort$rdis<100,"S","H")))
  Cohort$adl=ifelse(Cohort$radl <1,"A",ifelse(Cohort$radl<8,"L",ifelse(Cohort$radl<100,"L","H")))
  
  #### removing observation only have one dimension 
  
  Cohort$adl[which(is.na(Cohort$dis) & !is.na(Cohort$adl))] = NA
  Cohort$dis[which(is.na(Cohort$adl) & !is.na(Cohort$dis))] = NA
  
  #### combine two dimensions 
  Cohort$state = paste0(Cohort$dis,Cohort$adl)
  Cohort$state[which(Cohort$state=="NANA")] = NA
  unique(Cohort$state)
  Cohort$state = factor(Cohort$state,level= c("NA","NL","SA","SL","HH"))
  
  ### remove empty waves and fill age
  Cohort[which(Cohort$state=="HH"),"age"] = Cohort$radyear[which(Cohort$state=="HH")]-Cohort$rabyear[which(Cohort$state=="HH")]
  c = Cohort
  Mis = Cohort[which(is.na(Cohort$state)),]
  for (i in 1:nrow(Mis)) {
    I = Mis$rahhidpn[i]
    W = Mis$wave[i]
    if(length(is.na(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W+1),"state"]))==0|
       length(is.na(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W-1),"state"]))==0){next()}
    if(is.na(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W+1),"state"])|
       is.na(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W-1),"state"])){next()}
    if(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W+1),"state"]=="HH"){
      c[which(c$rahhidpn==I&c$wave==W),"state"]=Cohort[which(Cohort$rahhidpn==I&Cohort$wave==W-1),"state"]
      c[which(c$rahhidpn==I&c$wave==W),"age"] = ifelse((Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W-1)]+1)==(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1),(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1),
                                                       sample(c((Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W-1)]+1):(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1)),1))
    }else{
      n=rnorm(1)
      c[which(c$rahhidpn==I&c$wave==W),"state"]=ifelse(n>0,Cohort[which(Cohort$rahhidpn==I&Cohort$wave==W-1),"state"],Cohort[which(Cohort$rahhidpn==I&Cohort$wave==W+1),"state"])
      c[which(c$rahhidpn==I&c$wave==W),"age"] = ifelse((Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W-1)]+1)==(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1),(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1),
                                                       sample(c((Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W-1)]+1):(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1)),1))
    }
  }
  
  c_aw = left_join(c,dem)
  c_aw$pre_state = c_aw$state
  c_aw$state = lead(c_aw$state)
  c_aw$aw = ifelse(!is.na(c_aw$pre_state) & is.na(c_aw$state),1,0)
  c_aw$aw = ifelse(c_aw$pre_state=="HH",NA,c_aw$aw)
  c_aw$aw = ifelse(c_aw$wave==14,NA,c_aw$aw)
  denom.fit = glm(aw ~ as.factor(ragender) + as.factor(race) + as.factor(edu) + age + as.factor(pre_state), 
                  family = binomial(), data = c_aw, na.action = na.exclude)
  # pd <- predict(denom.fit, type = "response")
  
  numer.fit <- glm(aw~1, family = binomial(), data = c_aw, na.action = na.exclude)
  # pn <- predict(numer.fit, type = "response")
  
  Mis = c[which(is.na(c$state)),]
  c_hh=c()
  for (id in unique(c$rahhidpn)) {
    tem = c %>% filter(rahhidpn==id) 
    if(length(which(tem$state=="HH")[-1])>0){
      tem = tem[-which(tem$state=="HH")[-1],]
    }
    if(nrow(tem) >1){
      c_hh = rbind(c_hh,tem)
    }
  }
  
  ### tranform into state by age and fill gap ###
  c_con = c()
  for (id in unique(c_hh$rahhidpn)) {
    tem = c_hh %>% filter(rahhidpn==id)
    tem$wave = ifelse(tem$state=="HH",NA,tem$wave)
    if("TRUE" %in% duplicated(tem$age)){
      if(!"TRUE" %in% is.na(tem[which(tem$age == tem$age[which(duplicated(tem$age))]),"wave"]))
        tem = tem[-which(duplicated(tem$age)),]
    }
    tem = full_join(tem,data.frame(age= c(35:105)),by = "age")
    tem = arrange(tem,age)
    tem$rahhidpn=id
    Mis = tem[which(is.na(tem$state)),]
    for (i in 1:nrow(Mis)) {
      A = Mis$age[i]
      if(length(is.na(tem[which(tem$age==A-1),"state"]))==0|
         length(is.na(tem[which(tem$age==A+1),"state"]))==0|
         length(is.na(tem[which(tem$age==A+2),"state"]))==0){next()}
      if(is.na(tem[which(tem$age==A+1),"state"])[1] & 
         is.na(tem[which(tem$age==A+2),"state"])[1]){next()}
      if(is.na(tem[which(tem$age==A+1),"state"])[1] & 
         !is.na(tem[which(tem$age==A+2),"state"])[1]){
        tem[which(tem$age==A),"state"]=tem[which(tem$age==A-1),"state"]
      }else if("TRUE" %in% (tem[which(tem$age==A+1),"state"]=="HH")){
        tem[which(tem$age==A),"state"]=tem[which(tem$age==A-1),"state"]
      }else if(!is.na(tem[which(tem$age==A-1),"state"])&
               !is.na(tem[which(tem$age==A+1),"state"])){
        tem[which(tem$age==A),"state"]=sample(c(tem[which(tem$age==A-1),"state"],tem[which(tem$age==A+1),"state"]),1)
      }
    }
    c_con = rbind(c_con,tem)
  }
  c_con = c_con %>% filter(!is.na(age))
  
  c_con = c_con %>% group_by(rahhidpn) %>% mutate(pre_state = lag(state))
  c_con$dis = substring(c_con$state,0,1)
  c_con$adl = substring(c_con$state,2,3)
  c_con$pre_dis = substring(c_con$pre_state,0,1)
  c_con$pre_adl = substring(c_con$pre_state,2,3)
  # c_con = c_con %>% select(-raeduc,-ragender,-raracem,-edu,-rahispan,-race)
  
  
  c_con = left_join(c_con,dem)
  
  ### sample weight
  c_con = left_join(c_con,SW)
  c_con = c_con %>% group_by(rahhidpn) %>% mutate(lag_wt = lag(rwtcrnh)) %>% 
    mutate(wt = coalesce(rwtcrnh,lag_wt)) %>% mutate(lag_wt = lag(wt)) %>% 
    mutate(wt = coalesce(wt,lag_wt)) %>% mutate(lag_wt = lag(wt)) %>% 
    mutate(wt = coalesce(wt,lag_wt)) %>% select(-lag_wt,-rwtcrnh)
  
  c_con$pre_state = as.character(c_con$pre_state)
  c_con$pre_state = ifelse(c_con$pre_state =="HH",NA,c_con$pre_state)
  c_con$pre_state = factor(c_con$pre_state,levels = c("NA","NL","SA","SL"))
  pd <- predict(denom.fit, newdata = c_con, type = "response")
  pn <- predict(numer.fit, newdata = c_con, type = "response")
  c_con$w <- (1-pn)/(1-pd)
  
  
  save(list = c("Cohort","c_con"),file = paste0("output/C",substr(CO,3,4)," (10).Rdata"))
  
} 

# CMM----
## estimation ----
if(!dir.exists(file.path("output/CMM"))){dir.create(file.path("output/CMM"))}
YAW = data.frame(FY = c(1934,1944,1924,1934,1914,1924),
                 FA = c(60,60,70,70,80,80),
                 WA = c(4,9,4,9,4,9))

YAW$LY = YAW$FY+9
YAW$LA = YAW$FA+9  

popsize = 100000
# iteration = 1000
memory.limit(size = 1000000)
co = c("14","24","34","44")
load("output/dem.Rdata")
rm("Raw")
for (CO in co) {
  yaw = YAW %>% filter(FY==paste0("19",CO))
  for (row in 1:nrow(yaw)) {
    Age = yaw[row,"FA"]:(yaw[row,"FA"]+9)
    iteration = 100
    L = length(Age)
    for (z in 1:6) {
      if(z==6){iteration = 1}
      ALL = foreach(iter=1:iteration,.packages = c("tidyverse","nnet","haven")) %dopar% {
        
        if(iteration != 1){
        C_CON = c()
        for (COH in co) {
          load(paste0("output/C",COH," (10).RData"))
          tem = data.frame(rahhidpn = sample(unique(c_con$rahhidpn),length(unique(c_con$rahhidpn)),replace = T))
          c_con = left_join(tem,c_con)
          c_con$cohort = COH
          C_CON = rbind(C_CON,c_con)
        }
        }else{
          C_CON = c()
          for (COH in co) {
            load(paste0("output/C",COH," (10).RData"))
            c_con$cohort = COH
            C_CON = rbind(C_CON,c_con)
          }
        }
        
        BL = C_CON %>% filter(age%in% c((yaw[row,"FA"]-5):(yaw[row,"FA"]+4))) %>% filter(wave == yaw[row,"WA"])%>% filter(wt >0)
        
        BL = BL %>% group_by(state,ragender) %>% summarise(tot = sum(wt))
        BL = BL %>% filter(!is.na(ragender)&!is.na(state))
        BL = BL %>% mutate(state = as.character(state))
        BL = BL %>% ungroup() %>% mutate(pro = round(tot/sum(tot)*popsize,0)) %>% filter(pro>0)
        BL$state = recode(BL$state,"NA"=1,"NL"=2,"SA"=3,"SL"=4)
        
        bl = BL %>% mutate_at(c("ragender"),as.numeric)
        bl = as.matrix(bl[,c(1,2,4)])
        bl[,3] = as.numeric(bl[,3])
        
        
        tms = multinom(formula = state ~ pre_state + age + agesq + as.factor(ragender) + as.factor(cohort) +
                         as.factor(ragender):age + as.factor(cohort):age + as.factor(ragender):as.factor(cohort), 
                       weights = wt*w, data = C_CON %>% mutate(agesq = age^2),na.action=na.omit,maxit = 1000)
        
        dwrite <- crossing(pre_state = c("NA","NL","SA","SL"),cohort = co, ragender = c("1","2"),age=c(60:89))
        dwrite$agesq = dwrite$age ^ 2
        pp.writes <- cbind(dwrite, predict(tms, newdata = dwrite, type = "probs", se = TRUE))
        
        lpp.writes <- pp.writes %>% pivot_longer(c(6:10),names_to = "state",values_to = "prob")
        lpp.writes = lpp.writes %>% mutate(prob = ifelse(pre_state %in% c("SA","SL") & state %in% c("NA","NL"), 0, prob))
        lpp.writes$pre_state= recode(lpp.writes$pre_state,"NA"=1,"NL"=2,"SA"=3,"SL"=4)
        lpp.writes$state= recode(lpp.writes$state,"NA"=1,"NL"=2,"SA"=3,"SL"=4,"HH"=5)
        lpp.writes = lpp.writes %>% group_by(pre_state,ragender,age,cohort) %>% mutate(c = cumsum(prob))
        lpp.writes = lpp.writes %>% mutate(c = ifelse(state=="5", 1, c))
        lpp.writes = lpp.writes %>% select(-agesq,-prob)
        ## add in dead to dead prob
        HH = crossing(pre_state= 5,cohort=unique(lpp.writes$cohort),ragender= unique(lpp.writes$ragender),
                      age = unique(lpp.writes$age),state = c(1,2,3,4,5),c = 0)
        HH = HH %>% mutate(c = ifelse(pre_state %in% c(5) & state %in% c(5), 1, c))
        lpp.writes = bind_rows(lpp.writes,HH)
        lpp.writes = xtabs(c ~ pre_state+state+cohort+age+ragender,data = lpp.writes)
        
        # microsimulation
        Res = apply(bl, 1,function(x) matrix(c(rep(x[2], x[3]),rep(x[1], x[3])),nrow = x[3]))
        Res = do.call(rbind, Res)
        
        for(age in Age[-L]){
          Y = apply(Res, 1, function(x) lpp.writes[as.character(x[age-(Age[1]-2)]),,CO,as.character(age),as.character(x[1])])
          Y = t(Y)
          
          NEX = apply(Y,1,function(y){
            RAM = runif(1)
            s = ifelse(RAM<y[1],1,ifelse(RAM<y[2],2,ifelse(RAM<y[3],3,ifelse(RAM<y[4],4,5))))
            s
          })
          Res= cbind(Res,NEX)
        }
        Res
      }
      assign(paste0("ALL",CO,yaw[row,"FA"],".",z),ALL)
      save(list = c(paste0("ALL",CO,yaw[row,"FA"],".",z)),file = paste0("output/CMM/ALL",CO,yaw[row,"FA"],".",z,".RData"))
      rm(list = paste0("ALL",CO,yaw[row,"FA"],".",z))
    }
  }
}
stopImplicitCluster() 

co = c("14","24","34","44")
for (CO in co) {
  yaw = YAW %>% filter(FY==paste0("19",CO))
  for (row in 1:nrow(yaw)) {
    Age = yaw[row,"FA"]:(yaw[row,"FA"]+9)
    L = length(Age)
    for (z in 1:6) {
      load(paste0("output/CMM/ALL",CO,yaw[row,"FA"],".",z,".RData"))
      ALL=get(paste0("ALL",CO,yaw[row,"FA"],".",z))
      
      iteration = length(ALL)
      TIME<- vector(mode = "list", length = iteration)
      
      for(iter in 1:iteration){
        empty = crossing(id= c(1:nrow(ALL[[iter]])),state=c("NA", "NL", "SA", "SL", "HH"))
        empty$mob = substring(empty$state,0,1)
        empty$lim = substring(empty$state,2,3)
        
        S_all = as.data.frame(ALL[[iter]])
        colnames(S_all)=c("sex",Age,Age[L]+1)
        S_all <- S_all %>% 
          mutate_at(c(2:12), funs(recode(., `1`="NA", `2`="NL",`3`="SA",`4`="SL",`5`="HH")))
        S_all$id = 1:nrow(S_all)
        tem = S_all[,c(1,2,ncol(S_all))]
        colnames(tem)[2] = "ini"
        
        time_spend=pivot_longer(S_all,c(3:(2+L-1))) %>% count(id,value)
        time_spend2=pivot_longer(S_all,c(2)) %>% count(id,value) %>% mutate(n1=n/2) %>% select(-n)
        time_spend3=pivot_longer(S_all,c(2+L)) %>% count(id,value) %>% mutate(n2=n/2) %>% select(-n)
        
        time_spend = full_join(empty,time_spend,by = c("state"="value","id"))
        time_spend= full_join(time_spend,time_spend2,by = c("id", "state"="value"))
        time_spend= full_join(time_spend,time_spend3,by = c("id", "state"="value"))
        time_spend = time_spend %>% mutate(n= ifelse(is.na(n),0,n), n1=ifelse(is.na(n1),0,n1),n2=ifelse(is.na(n2),0,n2)) %>% 
          mutate(n=n+n1+n2) %>% select(-n1,-n2)
        
        time_spend = full_join(time_spend,tem,by = c("id"))
        
        
        # tem = time_spend %>% group_by(id,mob) %>% summarise(x = sum(n)) %>% filter(mob =="S") %>% mutate(df=ifelse(x==0,T,F))
        # time_spend = full_join(time_spend,tem[,c(1,4)],by = c("id"))
        # time_spend %>% group_by(sex,state) %>% summarise(Freq = mean(n))
        time_spend$iter = iter+(z-1)*100
        TIME[[iter]]=time_spend
      }
      assign(paste0("TIME",CO,yaw[row,"FA"],".",z),TIME)
      save(list = c(paste0("TIME",CO,yaw[row,"FA"],".",z)),file = paste0("output/CMM/TIME",CO,yaw[row,"FA"],".",z,".RData"))
      rm(list = c(paste0("TIME",CO,yaw[row,"FA"],".",z),"ALL","TIME","S_all","time_spend"))
    }
  }
}


## result ####
ver = "CMM"

library(tidyverse)

YAW = data.frame(FY = c(1934,1944,1924,1934,1914,1924),
                 FA = c(60,60,70,70,80,80),
                 WA = c(4,9,4,9,4,9))

YAW$LY = YAW$FY+9
YAW$LA = YAW$FA+9

SEX = function(df){df %>% group_by(mob,lim,sex,iter) %>% summarise(n = mean(n),.groups="drop")}
SEXin = function(df){df$ini = substring(df$ini,0,1)
df %>% group_by(mob,lim,sex,ini,iter) %>% summarise(n = mean(n),.groups="drop")}
RES = vector("list", nrow(YAW))
names(RES) <- paste0(substr(YAW[,"FY"],3,4),YAW[,"FA"])
for (a in unique(YAW$FA)) {
  yaw = YAW %>% filter(FA== a)
  for (row in 1:nrow(yaw)) {
    t = paste0("TIME",substr(yaw[row,"FY"],3,4),a,".")
    SEXdf = c()
    SEXidf = c()
    SEXindf = c()
    for(i in 1:6){
      load(paste0("output/",ver,"/",t,i,".RData"))
      tab = lapply(get(paste0(t,i)), SEX) %>% bind_rows()
      SEXdf = rbind(SEXdf,tab)
      tab = lapply(get(paste0(t,i)), SEXin) %>% bind_rows()
      SEXidf = rbind(SEXidf,tab)
      rm(list = c(paste0(t,i),"tab"))
    }
    res = lapply(ls(pattern="*df"), function(x) get(x))
    names(res) = ls(pattern="*df")
    RES[[paste0(substr(yaw[row,"FY"],3,4),a)]] = res
    rm(list = c(ls(pattern='*df'),"res"))
  }
}
save(list = c("RES"),file = paste0("output/",ver,"/RES.RData"))



# MMM (recursive) ----
## estimation ---- 
if(!dir.exists(file.path("output/MMM_rec"))){dir.create(file.path("output/MMM_rec"))}
YAW = data.frame(FY = c(1934,1944,1924,1934,1914,1924),
                 FA = c(60,60,70,70,80,80),
                 WA = c(4,9,4,9,4,9))

YAW$LY = YAW$FY+9
YAW$LA = YAW$FA+9  

popsize = 100000
# iteration = 1000
memory.limit(size = 1000000)
co = c("14","24","34","44")
load("output/dem.Rdata")
rm("Raw")
registerDoParallel(min(detectCores(),16))
for (CO in co) {
  yaw = YAW %>% filter(FY==paste0("19",CO))
  
  for (row in 1:nrow(yaw)) {
    Age = yaw[row,"FA"]:(yaw[row,"FA"]+9)
    iteration = 100
    L = length(Age)
    for (z in 1:6) {
      if(z==6){iteration = 1}
      ALL = foreach(iter=1:iteration,.packages = c("tidyverse","nnet","haven")) %dopar% {
        if(iteration != 1){
          C_CON = c()
          for (COH in co) {
            load(paste0("output/C",COH," (10).RData"))
            tem = data.frame(rahhidpn = sample(unique(c_con$rahhidpn),length(unique(c_con$rahhidpn)),replace = T))
            c_con = left_join(tem,c_con)
            c_con$cohort = COH
            C_CON = rbind(C_CON,c_con)
          }
        }else{
          C_CON = c()
          for (COH in co) {
            load(paste0("output/C",COH," (10).RData"))
            c_con$cohort = COH
            C_CON = rbind(C_CON,c_con)
          }
        }
        C_CON$pre_dis = factor(C_CON$pre_dis,levels = c("N","S"))
        C_CON$dis = factor(C_CON$dis,levels = c("N","S","H"))
        C_CON$pre_adl = factor(C_CON$pre_adl,levels = c("A","L"))
        C_CON$adl = factor(C_CON$adl,levels = c("A","L","H"))
        
        BL = C_CON %>% filter(age%in% c((yaw[row,"FA"]-5):(yaw[row,"FA"]+4))) %>% filter(wave == yaw[row,"WA"])%>% filter(wt >0)
        BL = BL %>% group_by(adl,dis,ragender) %>% summarise(tot = sum(wt))%>% filter(!is.na(ragender)&!is.na(adl)) %>% mutate(adl = as.character(adl))
        BL = BL %>% ungroup() %>% mutate(pro = round(tot/sum(tot)*popsize,0)) %>% filter(pro>0)
        
        BL$adl = recode(BL$adl,"A"=1,"L"=2)
        BL$dis = recode(BL$dis,"N"=1,"S"=2)
        
        bl = BL %>% mutate_at(c("ragender"),as.numeric)
        bl = as.matrix(bl[,c(1,2,3,5)])
        bl[,4] = as.numeric(bl[,4])
        
        
        tmsa = multinom(formula = adl ~ pre_adl+pre_dis+dis + age +pre_dis:age+pre_adl:age+dis:age+ agesq + as.factor(ragender)  + as.factor(cohort) +
                          as.factor(ragender):age + as.factor(ragender):pre_dis + as.factor(ragender):pre_adl+as.factor(ragender):dis+
                          as.factor(cohort):age + as.factor(ragender):as.factor(cohort),
                        weights = wt*w, data = C_CON %>% mutate(agesq = age^2) %>% filter(dis %in% c("N","S")),na.action=na.omit,maxit = 1000)
        
        dwrite <- crossing(pre_dis = c("N","S"),pre_adl = c("A","L"), dis = c("N","S"),cohort = co, ragender = c("1","2"),age=c(60:90))
        dwrite$agesq = dwrite$age ^ 2
        pp.writesa <- cbind(dwrite, predict(tmsa, newdata = dwrite, type = "probs", se = TRUE))
        colnames(pp.writesa)[ncol(pp.writesa)] = "L"
        pp.writesa$A = 1-pp.writesa$L
        pp.writesa$H = 0
        
        lpp.writesa <- pp.writesa %>% pivot_longer(c(8:10),names_to = "adl",values_to = "prob")
        lpp.writesa$pre_adl= recode(lpp.writesa$pre_adl,"A"=1,"L"=2)
        lpp.writesa$pre_dis= recode(lpp.writesa$pre_dis,"N"=1,"S"=2)
        lpp.writesa$dis = recode(lpp.writesa$dis,"N"=1,"S"=2)
        lpp.writesa$adl= recode(lpp.writesa$adl,"A"=1,"L"=2,"H"=3)
        lpp.writesa = lpp.writesa %>% group_by(dis,pre_dis,pre_adl,ragender,age,cohort) %>% arrange(adl) %>% mutate(c = cumsum(prob))
        lpp.writesa = lpp.writesa %>% mutate(c = ifelse(adl=="2", 1, c))
        lpp.writesa = lpp.writesa %>% select(-agesq,-prob)
        ## add in dead to dead prob
        H = crossing(pre_adl= c(1,2,3),pre_dis= c(1,2,3),dis = 3,cohort=unique(lpp.writesa$cohort),ragender= unique(lpp.writesa$ragender),
                     age = unique(lpp.writesa$age),adl = c(1,2,3),c = 0)
        H = H %>% mutate(c = ifelse(adl==3, 1, c))
        lpp.writesa = bind_rows(lpp.writesa,H)
        lpp.writesa = xtabs(c ~ pre_adl+adl+cohort+pre_dis+dis+age+ragender,data = lpp.writesa)
        
        
        tmsd = multinom(formula = dis ~ pre_dis +pre_adl + age +pre_adl:age+pre_dis:age+ agesq + as.factor(ragender) +as.factor(cohort) + 
                          as.factor(ragender):age + as.factor(cohort):age + as.factor(ragender):pre_dis + as.factor(ragender):pre_adl+
                          as.factor(ragender):as.factor(cohort),
                        weights = wt*w, data = C_CON %>% mutate(agesq = age^2) %>% filter(pre_adl %in%c("A","L")),na.action=na.omit,maxit = 1000)
        
        dwrite <- crossing(pre_dis = c("N","S"),pre_adl = c("A","L"),cohort = co, ragender = c("1","2"),age=c(60:90))
        dwrite$agesq = dwrite$age ^ 2
        pp.writesd <- cbind(dwrite, predict(tmsd, newdata = dwrite, type = "probs", se = TRUE))
        
        lpp.writesd <- pp.writesd %>% pivot_longer(c(7:9),names_to = "dis",values_to = "prob")
        lpp.writesd = lpp.writesd %>% mutate(prob = ifelse(pre_dis %in% c("S") & dis %in% c("N"), 0, prob))
        lpp.writesd$pre_adl= recode(lpp.writesd$pre_adl,"A"=1,"L"=2)
        lpp.writesd$pre_dis= recode(lpp.writesd$pre_dis,"N"=1,"S"=2)
        lpp.writesd$dis= recode(lpp.writesd$dis,"N"=1,"S"=2,"H"=3)
        lpp.writesd = lpp.writesd %>% group_by(pre_dis,pre_adl,ragender,age,cohort) %>% arrange(dis) %>% mutate(c = cumsum(prob))
        lpp.writesd = lpp.writesd %>% mutate(c = ifelse(dis=="3", 1, c))
        lpp.writesd = lpp.writesd %>% select(-agesq,-prob)
        ## add in dead to dead prob
        H = crossing(pre_adl= c(1,2,3),pre_dis = 3,cohort=unique(lpp.writesd$cohort),ragender= unique(lpp.writesd$ragender),
                     age = unique(lpp.writesd$age),dis = c(1,2,3),c = 0)
        H = H %>% mutate(c = ifelse(dis %in% c(3), 1, c))
        lpp.writesd = bind_rows(lpp.writesd,H)
        lpp.writesd = xtabs(c ~ pre_dis+dis+cohort+pre_adl+age+ragender,data = lpp.writesd)
        
        # microsimulation
        Resa = apply(bl, 1,function(x) matrix(c(rep(x[3], x[4]),rep(x[1], x[4])),nrow = x[4]))
        Resd = apply(bl, 1,function(x) matrix(c(rep(x[3], x[4]),rep(x[2], x[4])),nrow = x[4]))
        Resa = do.call(rbind, Resa)
        Resd = do.call(rbind, Resd)
        
        for(age in Age){
          resd = cbind(Resd,Resa[,ncol(Resa)])
          Yd = apply(resd, 1, function(x) lpp.writesd[as.character(x[age-(Age[1]-2)]),,CO,as.character(x[age-(Age[1]-3)]),as.character(age),as.character(x[1])])
          Yd = t(Yd)
          
          RAM = runif(nrow(Resd))
          Yd = cbind(Yd,RAM)
          NEXd = apply(Yd,1,function(y){
            s = ifelse(y[4]<y[1],1,ifelse(y[4]<y[2],2,3))
            s
          })
          Resd= cbind(Resd,NEXd)
          
          resa = cbind(Resa,Resd[,(ncol(Resd)-1):ncol(Resd)])
          Ya = apply(resa, 1, function(x) lpp.writesa[as.character(x[age-(Age[1]-2)]),,CO,as.character(x[age-(Age[1]-3)]),as.character(x[age-(Age[1]-4)]),as.character(age),as.character(x[1])])
          Ya = t(Ya)
          
          RAM = runif(nrow(Resa))
          Ya = cbind(Ya,RAM)
          NEXa = apply(Ya,1,function(y){
            s = ifelse(y[4]<y[1],1,ifelse(y[4]<y[2],2,3))
            s
          })
          Resa= cbind(Resa,NEXa)
          
        }
        
        Resa = as.data.frame(Resa)
        colnames(Resa)=c("sex",Age,Age[L]+1)
        Resa <- Resa %>% 
          mutate_at(c(2:12), funs(recode(., `1`="A", `2`="L",`3`="H")))
        Resa = as.matrix(Resa[,2:12])
        
        Resd = as.data.frame(Resd)
        colnames(Resd)=c("sex",Age,Age[L]+1)
        Resd <- Resd %>% 
          mutate_at(c(2:12), funs(recode(., `1`="N", `2`="S",`3`="H")))
        Res = Resd[,1]
        Resd = as.matrix(Resd[,2:12])
        
        tem = matrix(paste0(Resd,Resa),sum(bl[,4]),11)
        Res = cbind(Res,tem)
        
      }
      assign(paste0("ALL",CO,yaw[row,"FA"],".",z),ALL)
      save(list = c(paste0("ALL",CO,yaw[row,"FA"],".",z)),file = paste0("output/MMM_rec/ALL",CO,yaw[row,"FA"],".",z,".RData"))
      rm(list = paste0("ALL",CO,yaw[row,"FA"],".",z))
    }
  }
}
stopImplicitCluster() 

co = c("14","24","34","44")
for (CO in co) {
  yaw = YAW %>% filter(FY==paste0("19",CO))
  for (row in 1:nrow(yaw)) {
    Age = yaw[row,"FA"]:(yaw[row,"FA"]+9)
    L = length(Age)
    for (z in 1:6) {
      load(paste0("output/MMM_rec/ALL",CO,yaw[row,"FA"],".",z,".RData"))
      ALL=get(paste0("ALL",CO,yaw[row,"FA"],".",z))
      
      iteration = length(ALL)
      TIME<- vector(mode = "list", length = iteration)
      
      for(iter in 1:iteration){
        empty = crossing(id= c(1:nrow(ALL[[iter]])),state=c("NA", "NL", "SA", "SL", "HH"))
        empty$mob = substring(empty$state,0,1)
        empty$lim = substring(empty$state,2,3)
        
        S_all = as.data.frame(ALL[[iter]])
        colnames(S_all)=c("sex",Age,Age[L]+1)
        S_all <- S_all %>% 
          mutate_at(c(2:12), funs(recode(., `1`="NA", `2`="NL",`3`="SA",`4`="SL",`5`="HH")))
        S_all$id = 1:nrow(S_all)
        tem = S_all[,c(1,2,ncol(S_all))]
        colnames(tem)[2] = "ini"
        
        time_spend=pivot_longer(S_all,c(3:(2+L-1))) %>% count(id,value)
        time_spend2=pivot_longer(S_all,c(2)) %>% count(id,value) %>% mutate(n1=n/2) %>% select(-n)
        time_spend3=pivot_longer(S_all,c(2+L)) %>% count(id,value) %>% mutate(n2=n/2) %>% select(-n)
        
        time_spend = full_join(empty,time_spend,by = c("state"="value","id"))
        time_spend= full_join(time_spend,time_spend2,by = c("id", "state"="value"))
        time_spend= full_join(time_spend,time_spend3,by = c("id", "state"="value"))
        time_spend = time_spend %>% mutate(n= ifelse(is.na(n),0,n), n1=ifelse(is.na(n1),0,n1),n2=ifelse(is.na(n2),0,n2)) %>% 
          mutate(n=n+n1+n2) %>% select(-n1,-n2)
        
        time_spend = full_join(time_spend,tem,by = c("id"))
        
        
        # tem = time_spend %>% group_by(id,mob) %>% summarise(x = sum(n)) %>% filter(mob =="S") %>% mutate(df=ifelse(x==0,T,F))
        # time_spend = full_join(time_spend,tem[,c(1,4)],by = c("id"))
        # time_spend %>% group_by(sex,state) %>% summarise(Freq = mean(n))
        time_spend$iter = iter+(z-1)*100
        TIME[[iter]]=time_spend
      }
      assign(paste0("TIME",CO,yaw[row,"FA"],".",z),TIME)
      save(list = c(paste0("TIME",CO,yaw[row,"FA"],".",z)),file = paste0("output/MMM_rec/TIME",CO,yaw[row,"FA"],".",z,".RData"))
      rm(list = c(paste0("TIME",CO,yaw[row,"FA"],".",z),"ALL","TIME","S_all","time_spend"))
    }
  }
}

## result ####
ver = "MMM_rec"

library(tidyverse)

YAW = data.frame(FY = c(1934,1944,1924,1934,1914,1924),
                 FA = c(60,60,70,70,80,80),
                 WA = c(4,9,4,9,4,9))

YAW$LY = YAW$FY+9
YAW$LA = YAW$FA+9

SEX = function(df){df %>% group_by(mob,lim,sex,iter) %>% summarise(n = mean(n),.groups="drop")}
SEXin = function(df){df$ini = substring(df$ini,0,1)
df %>% group_by(mob,lim,sex,ini,iter) %>% summarise(n = mean(n),.groups="drop")}
RES = vector("list", nrow(YAW))
names(RES) <- paste0(substr(YAW[,"FY"],3,4),YAW[,"FA"])
for (a in unique(YAW$FA)) {
  yaw = YAW %>% filter(FA== a)
  for (row in 1:nrow(yaw)) {
    t = paste0("TIME",substr(yaw[row,"FY"],3,4),a,".")
    SEXdf = c()
    SEXidf = c()
    SEXindf = c()
    for(i in 1:6){
      load(paste0("output/",ver,"/",t,i,".RData"))
      tab = lapply(get(paste0(t,i)), SEX) %>% bind_rows()
      SEXdf = rbind(SEXdf,tab)
      tab = lapply(get(paste0(t,i)), SEXin) %>% bind_rows()
      SEXidf = rbind(SEXidf,tab)
      rm(list = c(paste0(t,i),"tab"))
    }
    res = lapply(ls(pattern="*df"), function(x) get(x))
    names(res) = ls(pattern="*df")
    RES[[paste0(substr(yaw[row,"FY"],3,4),a)]] = res
    rm(list = c(ls(pattern='*df'),"res"))
  }
}
save(list = c("RES"),file = paste0("output/",ver,"/RES.RData"))



# MMM (reduced-form) ----
## estimation ---- 
if(!dir.exists(file.path("output/MMM_red"))){dir.create(file.path("output/MMM_red"))}
YAW = data.frame(FY = c(1934,1944,1924,1934,1914,1924),
                 FA = c(60,60,70,70,80,80),
                 WA = c(4,9,4,9,4,9))

YAW$LY = YAW$FY+9
YAW$LA = YAW$FA+9  

popsize = 100000
# iteration = 1000
memory.limit(size = 1000000)
co = c("14","24","34","44")
load("output/dem.Rdata")
rm("Raw")
registerDoParallel(min(detectCores(),16))
for (CO in co) {
  yaw = YAW %>% filter(FY==paste0("19",CO))
  
  for (row in 1:nrow(yaw)) {
    Age = yaw[row,"FA"]:(yaw[row,"FA"]+9)
    iteration = 100
    L = length(Age)
    for (z in 1:6) {
      if(z==6){iteration = 1}
      ALL = foreach(iter=1:iteration,.packages = c("tidyverse","nnet","haven")) %dopar% {
        if(iteration != 1){
          C_CON = c()
          for (COH in co) {
            load(paste0("output/C",COH," (10).RData"))
            tem = data.frame(rahhidpn = sample(unique(c_con$rahhidpn),length(unique(c_con$rahhidpn)),replace = T))
            c_con = left_join(tem,c_con)
            c_con$cohort = COH
            C_CON = rbind(C_CON,c_con)
          }
        }else{
          C_CON = c()
          for (COH in co) {
            load(paste0("output/C",COH," (10).RData"))
            c_con$cohort = COH
            C_CON = rbind(C_CON,c_con)
          }
        }
        C_CON$pre_dis = factor(C_CON$pre_dis,levels = c("N","S"))
        C_CON$dis = factor(C_CON$dis,levels = c("N","S","H"))
        C_CON$pre_adl = factor(C_CON$pre_adl,levels = c("A","L"))
        C_CON$adl = factor(C_CON$adl,levels = c("A","L","H"))
        
        BL = C_CON %>% filter(age%in% c((yaw[row,"FA"]-5):(yaw[row,"FA"]+4))) %>% filter(wave == yaw[row,"WA"])%>% filter(wt >0)
        BL = BL %>% group_by(adl,dis,ragender) %>% summarise(tot = sum(wt))%>% filter(!is.na(ragender)&!is.na(adl)) %>% mutate(adl = as.character(adl))
        BL = BL %>% ungroup() %>% mutate(pro = round(tot/sum(tot)*popsize,0)) %>% filter(pro>0)
        
        BL$adl = recode(BL$adl,"A"=1,"L"=2)
        BL$dis = recode(BL$dis,"N"=1,"S"=2)
        
        bl = BL %>% mutate_at(c("ragender"),as.numeric)
        bl = as.matrix(bl[,c(1,2,3,5)])
        bl[,4] = as.numeric(bl[,4])
        
        
        
        tmsa = multinom(formula = adl ~ pre_adl+pre_dis + age+pre_dis:age +pre_adl:age+ agesq + as.factor(ragender) + as.factor(cohort) +
                          as.factor(ragender):age + as.factor(ragender):pre_dis + as.factor(ragender):pre_adl + as.factor(cohort):age + as.factor(ragender):as.factor(cohort),
                        weights = wt*w, data = C_CON %>% mutate(agesq = age^2),na.action=na.omit,maxit = 1000)
        
        dwrite <- crossing(pre_dis = c("N","S"),pre_adl = c("A","L"),cohort = co, ragender = c("1","2"),age=c(60:90))
        dwrite$agesq = dwrite$age ^ 2
        pp.writesa <- cbind(dwrite, predict(tmsa, newdata = dwrite, type = "probs", se = TRUE))
        
        
        lpp.writesa <- pp.writesa %>% pivot_longer(c(7:9),names_to = "adl",values_to = "prob")
        lpp.writesa$pre_adl= recode(lpp.writesa$pre_adl,"A"=1,"L"=2)
        lpp.writesa$pre_dis= recode(lpp.writesa$pre_dis,"N"=1,"S"=2)
        lpp.writesa$adl= recode(lpp.writesa$adl,"A"=1,"L"=2,"H"=3)
        lpp.writesa = lpp.writesa %>% group_by(pre_dis,pre_adl,ragender,age,cohort) %>% arrange(adl) %>% mutate(c = cumsum(prob))
        lpp.writesa = lpp.writesa %>% mutate(c = ifelse(adl=="3", 1, c))
        lpp.writesa = lpp.writesa %>% select(-agesq,-prob)
        ## add in dead to dead prob
        H = crossing(pre_dis= c(1,2),pre_adl = 3,cohort=unique(lpp.writesa$cohort),ragender= unique(lpp.writesa$ragender),
                     age = unique(lpp.writesa$age),adl = c(1,2,3),c = 0)
        H = H %>% mutate(c = ifelse(pre_adl %in% c(3) & adl %in% c(3), 1, c))
        lpp.writesa = bind_rows(lpp.writesa,H)
        lpp.writesa = xtabs(c ~ pre_adl+adl+cohort+pre_dis+age+ragender,data = lpp.writesa)
        
        
        tmsd = multinom(formula = dis ~ pre_adl+pre_dis+ age+pre_dis:age +pre_adl:age + agesq + as.factor(ragender) + as.factor(cohort) +
                          as.factor(ragender):age + as.factor(ragender):pre_dis + as.factor(ragender):pre_adl + as.factor(cohort):age + as.factor(ragender):as.factor(cohort),
                        weights = wt*w, data = C_CON %>% mutate(agesq = age^2) %>% filter(dis !="H"),na.action=na.omit,maxit = 1000)
        
        dwrite <- crossing(pre_dis = c("N","S"),pre_adl = c("A","L"),cohort = co, ragender = c("1","2"),age=c(60:90))
        dwrite$agesq = dwrite$age ^ 2
        pp.writesd <- cbind(dwrite, predict(tmsd, newdata = dwrite, type = "probs", se = TRUE))
        colnames(pp.writesd)[ncol(pp.writesd)] = "S"
        pp.writesd$S = ifelse(pp.writesd$pre_dis == "S", 1, pp.writesd$S)
        pp.writesd$N = 1-pp.writesd$S
        
        lpp.writesd <- pp.writesd %>% pivot_longer(c(7:8),names_to = "dis",values_to = "prob")
        lpp.writesd = lpp.writesd %>% mutate(prob = ifelse(pre_dis %in% c("S") & dis %in% c("N"), 0, prob))
        lpp.writesd$pre_adl= recode(lpp.writesd$pre_adl,"A"=1,"L"=2)
        lpp.writesd$pre_dis= recode(lpp.writesd$pre_dis,"N"=1,"S"=2)
        lpp.writesd$dis= recode(lpp.writesd$dis,"N"=1,"S"=2)
        lpp.writesd = lpp.writesd %>% group_by(pre_dis,pre_adl,ragender,age,cohort) %>% arrange(dis) %>% mutate(c = cumsum(prob))
        lpp.writesd = lpp.writesd %>% mutate(c = ifelse(dis=="2", 1, c))
        lpp.writesd = lpp.writesd %>% select(-agesq,-prob)
        ## add in dead to dead prob
        H = crossing(pre_dis= c(1,2),pre_adl = 3,cohort=unique(lpp.writesd$cohort),ragender= unique(lpp.writesd$ragender),
                     age = unique(lpp.writesd$age),dis = c(1,2),c = 0)
        H = H %>% mutate(c = ifelse(pre_dis %in% c(3), 1, c))
        lpp.writesd = bind_rows(lpp.writesd,H)
        lpp.writesd = xtabs(c ~ pre_dis+dis+cohort+pre_adl+age+ragender,data = lpp.writesd)
        
        # microsimulation
        Resa = apply(bl, 1,function(x) matrix(c(rep(x[3], x[4]),rep(x[1], x[4])),nrow = x[4]))
        Resd = apply(bl, 1,function(x) matrix(c(rep(x[3], x[4]),rep(x[2], x[4])),nrow = x[4]))
        Resa = do.call(rbind, Resa)
        Resd = do.call(rbind, Resd)
        
        for(age in Age){
          resa = cbind(Resa,Resd[,ncol(Resd)])
          Ya = apply(resa, 1, function(x) lpp.writesa[as.character(x[age-(Age[1]-2)]),,CO,as.character(x[age-(Age[1]-3)]),as.character(age),as.character(x[1])])
          Ya = t(Ya)
          
          resd = cbind(Resd,Resa[,ncol(Resa)])
          Yd = apply(resd, 1, function(x) lpp.writesd[as.character(x[age-(Age[1]-2)]),,CO,as.character(x[age-(Age[1]-3)]),as.character(age),as.character(x[1])])
          Yd = t(Yd)
          
          RAM = runif(nrow(Resa))
          Ya = cbind(Ya,RAM)
          NEXa = apply(Ya,1,function(y){
            s = ifelse(y[4]<y[1],1,ifelse(y[4]<y[2],2,3))
            s
          })
          Resa= cbind(Resa,NEXa)
          
          RAM = runif(nrow(Resd))
          Yd = cbind(Yd,RAM)
          NEXd = apply(Yd,1,function(y){
            s = ifelse(y[3]<y[1],1,2)
            s
          })
          Resd= cbind(Resd,NEXd)
        }
        Resd[which(Resa==3)]=3
        
        Resa = as.data.frame(Resa)
        colnames(Resa)=c("sex",Age,Age[L]+1)
        Resa <- Resa %>% 
          mutate_at(c(2:12), funs(recode(., `1`="A", `2`="L",`3`="H")))
        Resa = as.matrix(Resa[,2:12])
        
        Resd = as.data.frame(Resd)
        colnames(Resd)=c("sex",Age,Age[L]+1)
        Resd <- Resd %>% 
          mutate_at(c(2:12), funs(recode(., `1`="N", `2`="S",`3`="H")))
        Res = Resd[,1]
        Resd = as.matrix(Resd[,2:12])
        
        tem = matrix(paste0(Resd,Resa),sum(bl[,4]),11)
        Res = cbind(Res,tem)
        
      }
      assign(paste0("ALL",CO,yaw[row,"FA"],".",z),ALL)
      save(list = c(paste0("ALL",CO,yaw[row,"FA"],".",z)),file = paste0("output/MMM_red/ALL",CO,yaw[row,"FA"],".",z,".RData"))
      rm(list = paste0("ALL",CO,yaw[row,"FA"],".",z))
    }
  }
}
stopImplicitCluster() 

for (CO in co) {
  yaw = YAW %>% filter(FY==paste0("19",CO))
  for (row in 1:nrow(yaw)) {
    Age = yaw[row,"FA"]:(yaw[row,"FA"]+9)
    L = length(Age)
    for (z in 1:6) {
      load(paste0("output/MMM_red/ALL",CO,yaw[row,"FA"],".",z,".RData"))
      ALL=get(paste0("ALL",CO,yaw[row,"FA"],".",z))
      
      iteration = length(ALL)
      TIME<- vector(mode = "list", length = iteration)
      
      for(iter in 1:iteration){
        empty = crossing(id= c(1:nrow(ALL[[iter]])),state=c("NA", "NL", "SA", "SL", "HH"))
        empty$mob = substring(empty$state,0,1)
        empty$lim = substring(empty$state,2,3)
        
        S_all = as.data.frame(ALL[[iter]])
        colnames(S_all)=c("sex",Age,Age[L]+1)
        S_all <- S_all %>% 
          mutate_at(c(2:12), funs(recode(., `1`="NA", `2`="NL",`3`="SA",`4`="SL",`5`="HH")))
        S_all$id = 1:nrow(S_all)
        tem = S_all[,c(1,2,ncol(S_all))]
        colnames(tem)[2] = "ini"
        
        time_spend=pivot_longer(S_all,c(3:(2+L-1))) %>% count(id,value)
        time_spend2=pivot_longer(S_all,c(2)) %>% count(id,value) %>% mutate(n1=n/2) %>% select(-n)
        time_spend3=pivot_longer(S_all,c(2+L)) %>% count(id,value) %>% mutate(n2=n/2) %>% select(-n)
        
        time_spend = full_join(empty,time_spend,by = c("state"="value","id"))
        time_spend= full_join(time_spend,time_spend2,by = c("id", "state"="value"))
        time_spend= full_join(time_spend,time_spend3,by = c("id", "state"="value"))
        time_spend = time_spend %>% mutate(n= ifelse(is.na(n),0,n), n1=ifelse(is.na(n1),0,n1),n2=ifelse(is.na(n2),0,n2)) %>% 
          mutate(n=n+n1+n2) %>% select(-n1,-n2)
        
        time_spend = full_join(time_spend,tem,by = c("id"))
        
        
        # tem = time_spend %>% group_by(id,mob) %>% summarise(x = sum(n)) %>% filter(mob =="S") %>% mutate(df=ifelse(x==0,T,F))
        # time_spend = full_join(time_spend,tem[,c(1,4)],by = c("id"))
        # time_spend %>% group_by(sex,state) %>% summarise(Freq = mean(n))
        time_spend$iter = iter+(z-1)*100
        TIME[[iter]]=time_spend
      }
      assign(paste0("TIME",CO,yaw[row,"FA"],".",z),TIME)
      save(list = c(paste0("TIME",CO,yaw[row,"FA"],".",z)),file = paste0("output/MMM_red/TIME",CO,yaw[row,"FA"],".",z,".RData"))
      rm(list = c(paste0("TIME",CO,yaw[row,"FA"],".",z),"ALL","TIME","S_all","time_spend"))
    }
  }
}

## result ####
ver = "MMM_red"

library(tidyverse)

YAW = data.frame(FY = c(1934,1944,1924,1934,1914,1924),
                 FA = c(60,60,70,70,80,80),
                 WA = c(4,9,4,9,4,9))

YAW$LY = YAW$FY+9
YAW$LA = YAW$FA+9

SEX = function(df){df %>% group_by(mob,lim,sex,iter) %>% summarise(n = mean(n),.groups="drop")}
SEXin = function(df){df$ini = substring(df$ini,0,1)
df %>% group_by(mob,lim,sex,ini,iter) %>% summarise(n = mean(n),.groups="drop")}
RES = vector("list", nrow(YAW))
names(RES) <- paste0(substr(YAW[,"FY"],3,4),YAW[,"FA"])
for (a in unique(YAW$FA)) {
  yaw = YAW %>% filter(FA== a)
  for (row in 1:nrow(yaw)) {
    t = paste0("TIME",substr(yaw[row,"FY"],3,4),a,".")
    SEXdf = c()
    SEXidf = c()
    SEXindf = c()
    for(i in 1:6){
      load(paste0("output/",ver,"/",t,i,".RData"))
      tab = lapply(get(paste0(t,i)), SEX) %>% bind_rows()
      SEXdf = rbind(SEXdf,tab)
      tab = lapply(get(paste0(t,i)), SEXin) %>% bind_rows()
      SEXidf = rbind(SEXidf,tab)
      rm(list = c(paste0(t,i),"tab"))
    }
    res = lapply(ls(pattern="*df"), function(x) get(x))
    names(res) = ls(pattern="*df")
    RES[[paste0(substr(yaw[row,"FY"],3,4),a)]] = res
    rm(list = c(ls(pattern='*df'),"res"))
  }
}
save(list = c("RES"),file = paste0("output/",ver,"/RES.RData"))
