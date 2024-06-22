library(tidyverse)
library(nnet)
library(doParallel)
library(haven)
library(data.table)

# clean data----
# only need to run once, cleaned data will be stored in folder "output2"
# data.Rdata is the R data.frame converted from "randhrs1992_2018v1_STATA.zip"
# https://hrsdata.isr.umich.edu/data-products/rand-hrs-longitudinal-file-2018

load("HRS data.RData") # it is named as variable `Raw`


MA= Raw %>% select("rahhidpn",
                   paste0("r",c(2:14),"mstat")) 
MA = pivot_longer(MA,c(2:ncol(MA)))
MA$wave= parse_number(MA$name)
MA$name = gsub("[0-9]","",MA$name)
MA = pivot_wider(MA,names_from=name,values_from = value)
MA$rmstat = ifelse(MA$rmstat%in%c(1,2,3),0,ifelse(MA$rmstat%in%c(4,5,6),1,ifelse(MA$rmstat%in%c(7),2,NA)))

HE = Raw %>% select("rahhidpn",
                    paste0("r",c(2:14),"shlt")) 
HE = pivot_longer(HE,c(2:ncol(HE)))
HE$wave= parse_number(HE$name)
HE$name = gsub("[0-9]","",HE$name)
HE = pivot_wider(HE,names_from=name,values_from = value)

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
save.image(file = "output2/dem.Rdata")

for (WA in c(4,9)) {
  Period = Raw %>% select("rahhidpn","rabyear","radyear",paste0("r",c(WA:(WA+5)),"agey_e")) %>%
    mutate(WD = as.integer(radyear/2-995),WD = ifelse(is.na(WD),14,WD))
  Period = pivot_longer(Period,c(4:9),values_to = "age")
  Period$wave= parse_number(Period$name)
  Period = Period %>% select(-name)
  
  Period = left_join(Period,HE)
  Period$rshlt=ifelse(Period$wave>Period$WD,999,Period$rshlt)
  
  Period = left_join(Period,MA)
  Period$rmstat=ifelse(Period$wave>Period$WD,999,Period$rmstat)
  
  Period$mar=ifelse(Period$rmstat <1,"M",ifelse(Period$rmstat<2,"D",ifelse(Period$rmstat<100,"W","H")))
  Period$she=ifelse(Period$rshlt <3,"G",ifelse(Period$rshlt<4,"F",ifelse(Period$rshlt<6,"P","H")))
  
  #### removing obsevation only have one dimension 
  
  Period$mar[which(is.na(Period$she) & !is.na(Period$mar))] = NA
  Period$she[which(is.na(Period$mar) & !is.na(Period$she))] = NA
  
  #### combine two dimensions 
  Period$state = paste0(Period$mar,Period$she)
  Period$state[which(Period$state=="NANA")] = NA
  unique(Period$state)
  Period$state = factor(Period$state,level= c("MG","MF","MP","DG","DF","DP","WG","WF","WP","HH"))
  
  ### remove empty waves and fill age
  Period[which(Period$state=="HH"),"age"] = Period$radyear[which(Period$state=="HH")]-Period$rabyear[which(Period$state=="HH")]
  Period = Period %>% filter(age >=50)
  c = Period
  Mis = Period[which(is.na(Period$state)),]
  for (i in 1:nrow(Mis)) {
    I = Mis$rahhidpn[i]
    W = Mis$wave[i]
    if(length(is.na(Period[which(Period$rahhidpn == I & Period$wave==W+1),"state"]))==0|
       length(is.na(Period[which(Period$rahhidpn == I & Period$wave==W-1),"state"]))==0){next()}
    if(is.na(Period[which(Period$rahhidpn == I & Period$wave==W+1),"state"])|
       is.na(Period[which(Period$rahhidpn == I & Period$wave==W-1),"state"])){next()}
    if(Period[which(Period$rahhidpn == I & Period$wave==W+1),"state"]=="HH"){
      c[which(c$rahhidpn==I&c$wave==W),"state"]=Period[which(Period$rahhidpn==I&Period$wave==W-1),"state"]
      c[which(c$rahhidpn==I&c$wave==W),"age"] = ifelse((Period$age[which(Period$rahhidpn==I&Period$wave==W-1)]+1)==(Period$age[which(Period$rahhidpn==I&Period$wave==W+1)]-1),(Period$age[which(Period$rahhidpn==I&Period$wave==W+1)]-1),
                                                       sample(c((Period$age[which(Period$rahhidpn==I&Period$wave==W-1)]+1):(Period$age[which(Period$rahhidpn==I&Period$wave==W+1)]-1)),1))
    }else{
      n=rnorm(1)
      c[which(c$rahhidpn==I&c$wave==W),"state"]=ifelse(n>0,Period[which(Period$rahhidpn==I&Period$wave==W-1),"state"],Period[which(Period$rahhidpn==I&Period$wave==W+1),"state"])
      c[which(c$rahhidpn==I&c$wave==W),"age"] = ifelse((Period$age[which(Period$rahhidpn==I&Period$wave==W-1)]+1)==(Period$age[which(Period$rahhidpn==I&Period$wave==W+1)]-1),(Period$age[which(Period$rahhidpn==I&Period$wave==W+1)]-1),
                                                       sample(c((Period$age[which(Period$rahhidpn==I&Period$wave==W-1)]+1):(Period$age[which(Period$rahhidpn==I&Period$wave==W+1)]-1)),1))
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
    tem = full_join(tem,data.frame(age= c(50:105)),by = "age")
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
  c_con$mar = substring(c_con$state,0,1)
  c_con$she = substring(c_con$state,2,3)
  c_con$pre_mar = substring(c_con$pre_state,0,1)
  c_con$pre_she = substring(c_con$pre_state,2,3)
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
  c_con$pre_state = factor(c_con$pre_state,levels = c("MG","MF","MP","DG","DF","DP","WG","WF","WP"))
  pd <- predict(denom.fit, newdata = c_con, type = "response")
  pn <- predict(numer.fit, newdata = c_con, type = "response")
  c_con$w <- (1-pn)/(1-pd)
  
  save(list = c("Period","c_con"),file = paste0("output2/P",WA,"-",WA+5,".Rdata"))
  
} 


# estimation ---- 
popsize = 100000
# iteration = 1000
memory.limit(size = 1000000)
co = c(9)
a = 55
load("output2/dem.Rdata")
rm("Raw")
registerDoParallel(detectCores())
for (CO in co) {
  iteration = 100
  load(paste0("output2/P",CO,"-",CO+5,".Rdata"))
  BL = c_con %>% filter(age %in% c((a-4):(a+5)),wt >0) %>% group_by(rahhidpn) %>% arrange(age) %>% slice_head(n=1)
  for (z in 1:5) {
    ALL = foreach(iter=1:iteration,.packages = c("tidyverse","nnet","haven")) %dopar% {
      
      if(iteration==1){
        tem = data.frame(rahhidpn = unique(c_con$rahhidpn))
      }else{
        tem = data.frame(rahhidpn = sample(unique(c_con$rahhidpn),length(unique(c_con$rahhidpn)),replace = T))
      }
      C_CON = left_join(tem,c_con)
      
      C_CON$pre_mar = factor(C_CON$pre_mar,levels = c("M","D","W"))
      C_CON$mar = factor(C_CON$mar,levels = c("M","D","W","H"))
      C_CON$pre_she = factor(C_CON$pre_she,levels = c("G","F","P"))
      C_CON$she = factor(C_CON$she,levels = c("G","F","P","H"))
      C_CON = C_CON %>% filter(!(pre_mar == "D" & mar == "W"))
      C_CON = C_CON %>% filter(!(pre_mar == "W" & mar == "D"))
      
      Baseline = left_join(tem,BL) %>% filter(!is.na(wt),!mar == "H")
      Baseline = Baseline %>% group_by(she,mar,ragender) %>% summarise(tot = sum(wt))%>% filter(!is.na(ragender)&!is.na(she)) %>% mutate(she = as.character(she))
      Baseline = Baseline %>% ungroup() %>% mutate(pro = round(tot/sum(tot)*popsize,0)) %>% filter(pro>0)
      
      Baseline$she = recode(Baseline$she,"G"=1,"F"=2,"P"=3)
      Baseline$mar = recode(Baseline$mar,"M"=1,"D"=2,"W"=3)
      
      bl = Baseline %>% mutate_at(c("ragender"),as.numeric)
      bl = as.matrix(bl[,c(1,2,3,5)])
      bl[,4] = as.numeric(bl[,4])
      
      
      tmsa = multinom(formula = she ~ pre_she + pre_mar + age + agesq + as.factor(ragender) + as.factor(ragender):age,
                      weights = wt*w, data = C_CON %>% mutate(agesq = age^2) %>% filter(pre_she!="H"),na.action=na.omit,maxit = 1000)
      
      dwrite <- crossing(pre_mar = c("M","D","W"),pre_she = c("G","F","P"), ragender = c("1","2"),age=c(55:105))
      dwrite$agesq = dwrite$age ^ 2
      pp.writesa <- cbind(dwrite, predict(tmsa, newdata = dwrite, type = "probs", se = TRUE))
      
      lpp.writesa <- pp.writesa %>% pivot_longer(c(6:9),names_to = "she",values_to = "prob")
      lpp.writesa$pre_she= recode(lpp.writesa$pre_she,"G"=1,"F"=2,"P"=3)
      lpp.writesa$pre_mar= recode(lpp.writesa$pre_mar,"M"=1,"D"=2,"W"=3)
      lpp.writesa$she= recode(lpp.writesa$she,"G"=1,"F"=2,"P"=3,"H"=4)
      lpp.writesa = lpp.writesa %>% group_by(pre_mar,pre_she,ragender,age) %>% mutate(prob = 1/sum(prob)*prob)
      lpp.writesa = lpp.writesa %>% group_by(pre_mar,pre_she,ragender,age) %>% arrange(she) %>% mutate(c = cumsum(prob))
      lpp.writesa = lpp.writesa %>% mutate(c = ifelse(she=="4", 1, c))
      lpp.writesa = lpp.writesa %>% select(-agesq,-prob)
      ## add in dead to dead prob
      H = crossing(pre_mar= c(1,2,3),pre_she = 4,ragender= unique(lpp.writesa$ragender),
                   age = unique(lpp.writesa$age),she = c(1,2,3,4),c = 0)
      H = H %>% mutate(c = ifelse(pre_she %in% c(4) & she %in% c(4), 1, c))
      lpp.writesa = bind_rows(lpp.writesa,H)
      lpp.writesa = xtabs(c ~ pre_she+she+pre_mar+age+ragender,data = lpp.writesa)
      
      
      tmsd = multinom(formula = mar ~ pre_mar + age + agesq + as.factor(ragender) + as.factor(ragender):age ,
                      weights = wt*w, data = C_CON %>% mutate(agesq = age^2) %>% filter(mar !="H"),na.action=na.omit,maxit = 1000)
      
      dwrite <- crossing(pre_mar = c("M","D","W"), ragender = c("1","2"),age=c(55:105))
      dwrite$agesq = dwrite$age ^ 2
      pp.writesd <- cbind(dwrite, predict(tmsd, newdata = dwrite, type = "probs", se = TRUE))
      
      lpp.writesd <- pp.writesd %>% pivot_longer(c(5:7),names_to = "mar",values_to = "prob")
      lpp.writesd$pre_mar= recode(lpp.writesd$pre_mar,"M"=1,"D"=2,"W"=3)
      lpp.writesd$mar= recode(lpp.writesd$mar,"M"=1,"D"=2,"W"=3)
      lpp.writesd = lpp.writesd %>% mutate(prob = ifelse(pre_mar=="2"&mar=="3", 0, prob)) %>% mutate(prob = ifelse(pre_mar=="3"&mar=="2", 0, prob))
      lpp.writesd = lpp.writesd %>% group_by(pre_mar,ragender,age) %>% mutate(prob = 1/sum(prob)*prob)
      lpp.writesd = lpp.writesd %>% group_by(pre_mar,ragender,age) %>% arrange(mar) %>% mutate(c = cumsum(prob))

      lpp.writesd = lpp.writesd %>% select(-agesq,-prob)

      lpp.writesd = xtabs(c ~ pre_mar+mar+age+ragender,data = lpp.writesd)
      
      # microsimulation
      Resa = apply(bl, 1,function(x) matrix(c(rep(x[3], x[4]),rep(x[1], x[4])),nrow = x[4]))
      Resd = apply(bl, 1,function(x) matrix(c(rep(x[3], x[4]),rep(x[2], x[4])),nrow = x[4]))
      Resa = do.call(rbind, Resa)
      Resd = do.call(rbind, Resd)
      
      for(age in a:104){
        resa = cbind(Resa,Resd[,ncol(Resd)])
        Ya = apply(resa, 1, function(x) lpp.writesa[as.character(x[age-(a-2)]),,as.character(x[age-(a-3)]),as.character(age),as.character(x[1])])
        Ya = t(Ya)
        
        resd = cbind(Resd,Resa[,ncol(Resa)])
        Yd = apply(resd, 1, function(x) lpp.writesd[as.character(x[age-(a-2)]),,as.character(age),as.character(x[1])])
        Yd = t(Yd)
        
        RAM = runif(nrow(Resa))
        Ya = cbind(Ya,RAM)
        NEXa = apply(Ya,1,function(y){
          # change y[z], z is the number of categories + 1
          s = ifelse(y[5]<y[1],1,ifelse(y[5]<y[2],2,ifelse(y[5]<y[3],3,4)))
          s
        })
        Resa= cbind(Resa,NEXa)
        
        RAM = runif(nrow(Resd))
        Yd = cbind(Yd,RAM)
        NEXd = apply(Yd,1,function(y){
          s = ifelse(y[4]<y[1],1,ifelse(y[4]<(y[2]),2,3))
          s
        })
        Resd= cbind(Resd,NEXd)
      }
      # sync death to the other vector
      Resd[which(Resa==4)]=4
      
      Resa = as.data.frame(Resa)
      colnames(Resa)=c("sex",a:105)
      Resa <- Resa %>% 
        mutate_at(c(2:ncol(Resa)), funs(recode(., `1`="G", `2`="F",`3`="P",`4`="H")))
      Resa = as.matrix(Resa[,2:ncol(Resa)])
      
      Resd = as.data.frame(Resd)
      colnames(Resd)=c("sex",a:105)
      Resd <- Resd %>% 
        mutate_at(c(2:ncol(Resd)), funs(recode(., `1`="M", `2`="D",`3`="W",`4`="H")))
      Res = Resd[,1]
      Resd = as.matrix(Resd[,2:ncol(Resd)])
      
      tem = matrix(paste0(Resd,Resa),sum(bl[,4]),ncol(Resd))
      Res = cbind(Res,tem)
      
    }
    assign(paste0("ALL",CO,".",z),ALL)
    save(list = c(paste0("ALL",CO,".",z)),file = paste0("output2/ALL",CO,".",z,".RData"))
    rm(list = paste0("ALL",CO,".",z))
  }
}


# results ----
res = c()
for (co in c(9)) {
  for (z in 1:5) {
    load(paste0("output2/ALL",co,".",z,".RData"))
    for (iter in 1:100) {
      
      tem = get(paste0("ALL",co,".",z))[[iter]]
      tem = as.data.frame(tem)
      tem$ini = tem[,2]
      colnames(tem) = c("Sex",55:105,"INI")
      tem$id = 1:nrow(tem)
      
      tem = tem %>% pivot_longer(2:52,names_to="Age")
      tem$mar = substring(tem$value,0,1)
      tem$she = substring(tem$value,2,3)
      tem$ini = substring(tem$INI,0,1)
      tem = setDT(tem)
      
      
      t1 = 55
      t = 65
      list_m64 = tem[Age %in% c(t1:t),.N,by=.(id,Sex, mar,ini)][mar=="M"&N==(t-t1+1),id]
      list_M64 = tem[Age %in% c(t1:t),.N,by=.(id,Sex, mar,ini)][ini=="M",id]
      list_H = tem[id %in% list_M64 & Age %in% c((t+1):105),.N,by=.(id,Sex, mar)][mar=="H"&N==(105-t),id]
      list_M64 = setdiff(list_M64,list_H)
      list_M64 = setdiff(list_M64,list_m64)
      tem2 = copy(tem)
      tem2[,freq := 1][Age==t|Age==105,freq := 0.5]
      
      empty = crossing(id = list_M64, she = c("G","F","P","H")) %>% left_join(tem[id %in% list_M64 & Age %in% c(t+1),.(id,Sex)])
      MX = tem2[id%in% list_M64 & Age %in% c(t:105),.(N=sum(freq)),by=.(id,she)]
      MX = left_join(empty,MX) %>% setDT()
      MX$N = fifelse(is.na(MX$N),0,MX$N)
      MX = MX[,.(chg = mean(N)),by=.(Sex,she)]
      
      empty = crossing(id = list_m64, she = c("G","F","P","H")) %>% left_join(tem[id %in% list_m64 & Age %in% c(t+1),.(id,Sex)])
      MM = tem2[id%in% list_m64 & Age %in% c(t:105),.(N=sum(freq)),by=.(id,she)]
      MM = left_join(empty,MM) %>% setDT()
      MM$N = fifelse(is.na(MM$N),0,MM$N)
      MM = MM[,.(con = mean(N)),by=.(Sex,she)]
      
      data = left_join(MM,MX)
      data$iter = iter+100*(z-1)
      data$p = co
      data$age = paste0(t1,"-",t)
      res = rbind(res,data)
      
      t1 = 65
      t = 75
      list_m74 = tem[id %in% list_m64 & Age %in% c(t1:t),.N,by=.(id,Sex, mar,ini)][mar=="M" & N==(t-t1+1),id]
      list_M74 = setdiff(list_m64,list_m74)
      list_H = tem[id %in% list_M74 & Age %in% c((t+1):105),.N,by=.(id,Sex, mar)][mar=="H"& N == (105-t),id]
      list_M74 = setdiff(list_M74,list_H)
      tem2 = copy(tem)
      tem2[,freq := 1][Age==t|Age==105,freq := 0.5]
      
      empty = crossing(id = list_M74, she = c("G","F","P","H")) %>% left_join(tem[id %in% list_M74 & Age %in% c(t+1),.(id,Sex)])
      MX = tem2[id%in% list_M74 & Age %in% c(t:105),.(N=sum(freq)),by=.(id,she)]
      MX = left_join(empty,MX) %>% setDT()
      MX$N = fifelse(is.na(MX$N),0,MX$N)
      MX = MX[,.(chg = mean(N)),by=.(Sex,she)]
      
      empty = crossing(id = list_m74, she = c("G","F","P","H")) %>% left_join(tem[id %in% list_m74 & Age %in% c(t+1),.(id,Sex)])
      MM = tem2[id%in% list_m74 & Age %in% c(t:105),.(N=sum(freq)),by=.(id,she)]
      MM = left_join(empty,MM) %>% setDT()
      MM$N = fifelse(is.na(MM$N),0,MM$N)
      MM = MM[,.(con = mean(N)),by=.(Sex,she)]
      
      data = left_join(MM,MX)
      data$iter = iter+100*(z-1)
      data$p = co
      data$age = paste0(t1,"-",t)
      res = rbind(res,data)
      
      t1 = 75
      t = 85
      list_m84 = tem[id %in% list_m74 & Age %in% c(t1:t),.N,by=.(id,Sex, mar,ini)][mar=="M" & N==(t-t1+1),id]
      list_M84 = setdiff(list_m74,list_m84)
      list_H = tem[id %in% list_M84 & Age %in% c((t+1):105),.N,by=.(id,Sex, mar)][mar=="H"& N == (105-t),id]
      list_M84 = setdiff(list_M84,list_H)
      tem2 = copy(tem)
      tem2[,freq := 1][Age==t|Age==105,freq := 0.5]
      
      empty = crossing(id = list_M84, she = c("G","F","P","H")) %>% left_join(tem[id %in% list_M84 & Age %in% c(t+1),.(id,Sex)])
      MX = tem2[id%in% list_M84 & Age %in% c(t:105),.(N=sum(freq)),by=.(id,she)]
      MX = left_join(empty,MX) %>% setDT()
      MX$N = fifelse(is.na(MX$N),0,MX$N)
      MX = MX[,.(chg = mean(N)),by=.(Sex,she)]
      
      empty = crossing(id = list_m84, she = c("G","F","P","H")) %>% left_join(tem[id %in% list_m84 & Age %in% c(t+1),.(id,Sex)])
      MM = tem2[id%in% list_m84 & Age %in% c(t:105),.(N=sum(freq)),by=.(id,she)]
      MM = left_join(empty,MM) %>% setDT()
      MM$N = fifelse(is.na(MM$N),0,MM$N)
      MM = MM[,.(con = mean(N)),by=.(Sex,she)]
      
      data = left_join(MM,MX)
      data$iter = iter+100*(z-1)
      data$p = co
      data$age = paste0(t1,"-",t)
      res = rbind(res,data)
      
    }
    rm(list = c(paste0("ALL",co,".",z)))
  }
}
write_csv(res,paste0("output2/res.csv"))

Res = res %>% group_by(Sex,she,age,p) %>% summarise(m_con = median(con),l_con = quantile(con,0.025),u_con = quantile(con,0.975),
                                                    m_chg = median(chg),l_chg = quantile(chg,0.025),u_chg = quantile(chg,0.975))
write_csv(Res,paste0("output2/Summary.csv"))

