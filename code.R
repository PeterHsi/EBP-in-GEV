# package load ----
# install.packages(c("extRemes","Rfast"))
library("extRemes") # main package of EVA
library("Rfast")
library("boot")
library("matrixStats")
library("ggplot2")
library("tidyverse")
library("latex2exp")
library("readr")
library("see")
options(digits = 6)

# Some utility ------------------------------------------------------------
# 讀取資料夾內檔案名稱
gen_rrfnl <- function(cd, cd_layer=1){
    row_result_fn <- cd
    for(layer in 1:cd_layer){
        temp_rrfn <- NULL
        for(file_num in 1:length(row_result_fn)){
            temp_rrfn <- append(temp_rrfn,
                                paste(row_result_fn[file_num],
                                      list.files(row_result_fn[file_num]),
                                      sep="/"))
        }
        row_result_fn <- temp_rrfn
    }
    return(row_result_fn)
}

# 將 extReme::fevd 之結果(obj)讀取出來
ext_res <- function(obj, nd_s=F){
    par_mu <- obj$results$par
    NLL <- obj$results$value
    
    if (obj$par.models$location == "~t0" | obj$par.models$location == "~t+t0"){
        k <- length(obj$results$par)+1
    } else {
        k <- length(obj$results$par)
    }
    AIC <- 2*NLL+2*k
    BIC <- 2*NLL+log(length(obj$x))*k
    if (nd_s){
        par_sigma <- sqrt(diag(solve(obj$results$hessian)))
        list(par_mu = par_mu, par_sigma = par_sigma,
             NLL = NLL, AIC = AIC, BIC = BIC)
    } else {
        list(par_mu = par_mu, NLL = NLL, AIC = AIC, BIC = BIC)
    }
}

# Gen. GEV dist RV with breakpoint ----------------------------------------
# 論文假設的位置參數函數 mu_fun()
#     t 為斷點所在空間，a1 為參數函數對 t 之一次項係數， a2 為斷點處參數函數
#     斷層高度，tt 為斷點場所
mu_fun <- function(t, a0 = 0, a1 = 0, a2 = 0, tt = 0){a0+a1*t+a2*(t/max(t)>=tt)}

# 產生設定分配之隨機變數 tg_revd()
#     a1, a2, tt 沿用 mu_fun()，xi 為形狀參數設定，n 為樣本數
#     t 會自動產生並常規化至 [0, 1] 區間
tg_revd <- function(a1, a2, tt, xi, n){
    t <- c(1:n)/n
    mu <- mu_fun(t, a1 = a1, a2 = a2, tt = tt)
    x <- revd(n, loc = mu, scale = 1, shape = xi,
              type = "GEV")
    data <- data.frame(x, t)
    return(data)
}

# Breakpoint GEV MLE estimate ---------------------------------------------
# 帶有斷點的 GEV 估計(給定單一參數函數設定) fevd_bpt()
#     給定參數函數尋找斷點，data 輸入需為 dataframe
#     斷點所位於之變數必須命名為 t ，極端值需命名為 x
#     locf, scaf 為位置函數與規模函數，斷點須以 t0 命名
fevd_bpt <- function(data, locf="~t+t0", scaf="~1") {
    Break <- sort(unique(data$t))
    Break <- Break[2:(length(Break))]
    d <- numeric(length(Break))
    for (i in 1:length(Break)) {
        data$t0 <- (data$t >= Break[i])*1
        fit_temp <- fevd(x, data, location.fun=as.formula(locf),
                         scale.fun=as.formula(scaf), method="MLE", type="GEV")
        d[i] = fit_temp$results$value
    }
    bp <- Break[which.min(d)]
    return(bp)
}

# 報告帶有斷點的 GEV 估計(給定一群參數函數設定) fit_result()
#     data 之規定沿用 fevd_bpt()，locf_set, scaf_set 可為多個不同參數函數
#     程式會把所有組合估計完後，讀取估計最小負對數概似值，以此計算 AIC、BIC
#     並記錄各模型估計斷點位置
fit_result <- function(data, locf_set, scaf_set ="~1"){
    I <- length(locf_set); J <- length(scaf_set)
    NLL <- rep(0, I*J); AIC <- NLL; BIC <- NLL; set <- NLL; tc_est <- NLL
    for (j in 1:J){
        for (i in 1:I){
            if (locf_set[i] == "~t0" | locf_set[i] == "~t+t0"){
                tc <- fevd_bpt(data, locf_set[i])
                data$t0 <- (data$t >= tc)*1
            } else {
                tc <- NA
            }
            fit_temp <- fevd(x, data,
                             location.fun = as.formula(locf_set[i]),
                             scale.fun = as.formula(scaf_set[j]))
            set[(j-1)*I+i] <- paste(c(locf_set[i],scaf_set[j]), collapse="/")
            list_temp <- ext_res(fit_temp)
            NLL[(j-1)*I+i] <- list_temp$NLL
            AIC[(j-1)*I+i] <- list_temp$AIC
            BIC[(j-1)*I+i] <- list_temp$BIC
            tc_est[(j-1)*I+i] <- tc
        }
    }
    idx <- data.frame(set, NLL, AIC, BIC, tc_est)
    return(idx)
}

# Main --------------------------------------------------------------------

# 模擬函數 simulation()
#     一個參數設定會產生一個 csv 檔，檔頭加綴為 csvname
#     內容為各次迭代之斷點場所估計以及斷點對所有候選模型之 NLL、AIC、BIC 估計
#     a1_set, a2_set, tt_set, xi_set, n_set 為測試參數組合
#         (項量化 tg_revd() 之 a1, a2, tt, xi, n)
#     N 為迭代次數
#     locf_set, scaf_set 為估計之候選模型之參數函數
#         (沿用 fit_result()之 locf_set, scaf_set )
simulation <- function(a1_set, a2_set, tt_set, xi_set,
                       n_set, N, locf_set, scaf_set, csvname="result"){
    I <- length(locf_set); J <- length(scaf_set); K <- length(n_set)
    sl_NLL <- rep(0, N); sl_AIC <- sl_NLL; sl_BIC <- sl_NLL
    for (type4 in 1:length(xi_set)){
        for (type3 in 1:length(tt_set)){
            for (type2 in 1:length(a2_set)){
                for (type1 in 1:length(a1_set)){
                    for (k in 1:K) {
                        n <- n_set[k]; a1 <- a1_set[type1]; a2 <- a2_set[type2]
                        tt <- tt_set[type3]; xi <- xi_set[type4]
                        fname <- paste(gsub("[.]","p",paste(csvname,a1, a2, tt,
                                                            xi, n, sep="_")),
                                       ".csv",sep="")
                        first_counter <- TRUE
                        if ((tt == 0 & a2 != 0)|(a2 == 0 & tt != 0)){next}
                        for (i in 1:N){
                            data <- tg_revd(a1, a2, tt, xi, n)
                            fitt <- fit_result(data, locf_set, scaf_set)
                            if (first_counter){
                                header <- t(c("a1", "a2", "tt", "xi", "n", "i",
                                              paste("NLL:", fitt$set),
                                              paste("tc_est:", fitt$set)))
                                write.table(header, file = fname, sep = ", ",
                                            append = T, col.names = F, row.names = F)
                                first_counter = F
                            }
                            result <- t(c(a1,a2,tt,xi,n,i,fitt$NLL,fitt$tc_est))
                            write.table(result, file = fname, sep = ", ",
                                        append = T, col.names = F, row.names = F)
                        }
                    }
                }
            }
        }
    }
}

# 報告函數 report()
#     將模擬函數 simulation() 產生之 csv 檔產生報告之函數
#     cd, cd_layer 為指定檔案位置，會將路徑 cd 之下面第 cd_layer 層之資料夾內
#         所有 csv 讀取並彙整其：
#         - 給定參數設定。
#         - N 次估計各模型選擇準則選擇各模型次數
#         - 各模型估計斷點場所平均值與變異數
#     par_pos 此 list 為各設定參數在 csv 的 colmun
#     t0_pos 此 list 為各候選模型估計斷點場所在 csv 的 colmun
#     NLL_pos 此 list 為各候選模型之 NLL 估計值在 csv 的 colmun
#     n_pos 為樣本數在 csv 的 colmun
#     k 此 list 為各候選模型之使用的參數個數
#     tuedep 為利用給定參數(par_pos 指定讀取內容之前兩列是否為零)自動標籤給定模型
#     filenm 為處存之報告 csv 檔名
report <- function(cd, par_pos, t0_pos, NLL_pos, n_pos, k, tuedep,
                   cd_layer=1, filenm = "result"){
    # Get t0 pos info
    t0nam <- colnames(t0_pos)
    t0num <- as.array(as.matrix(t0_pos)[1:length(t0_pos)])
    NLLnam <- colnames(NLL_pos)
    NLLnum <- as.array(as.matrix(NLL_pos)[1:length(NLL_pos)])
    parnam <- colnames(par_pos)
    parnum <- as.array(as.matrix(par_pos)[1:length(par_pos)])
    
    # Get file info
    rrfnl <- gen_rrfnl(cd, cd_layer)
    
    # Loop
    total_result <- NULL
    for(i in 1:length(rrfnl)){
        temp <- as.data.frame(read.csv(rrfnl[i], row.names=NULL))
        temp_t0m <- as.vector(sapply(temp[t0num], mean, na.rm = T))
        temp_t0v <- as.vector(sapply(temp[t0num], var, na.rm = T))
        names(temp_t0m)=paste("t0m_", t0nam, sep="")
        names(temp_t0v)=paste("t0v_", t0nam, sep="")
        
        temp_par <- as.matrix(temp[1,parnum])
        names(temp_par)<-parnam
        whmol <- (temp_par==0)[1:2]==tuedep
        whmol <- which(whmol[1,] & whmol[2,])
        
        temp_NLL <- as.matrix(temp[NLLnum])
        colnames(temp_NLL)<-NLLnam
        temp_AIC <- t(2*t(temp_NLL)+2*k)
        temp_BIC <- t(2*t(temp_NLL)+log(temp[n_pos][[1,1]])*k)
        
        temp_NLLsl <- colSums(rowMins(temp_NLL)==temp_NLL)
        temp_AICsl <- colSums(rowMins(temp_AIC)==temp_AIC)
        temp_BICsl <- colSums(rowMins(temp_BIC)==temp_BIC)
        names(temp_NLLsl)<-paste("NLL",NLLnam,sep="_")
        names(temp_AICsl)<-paste("AIC",NLLnam,sep="_")
        names(temp_BICsl)<-paste("BIC",NLLnam,sep="_")
        
        tt0 <- switch(whmol,
                      t0mod <- NA,
                      t0mod <- 1,
                      t0mod <- NA,
                      t0mod <- 2)
        
        if(is.na(t0mod)){
            turt0m <- NA
            turt0v <- NA
        } else {
            turt0m <- temp_t0m[t0mod]
            turt0v <- temp_t0v[t0mod]
        }
        
        turMod <- names(tuedep[whmol])
        turNLLslr <- temp_NLLsl[whmol]/sum(temp_NLLsl)
        turAICslr <- temp_AICsl[whmol]/sum(temp_AICsl)
        turBICslr <- temp_BICsl[whmol]/sum(temp_BICsl)
        names(turt0m) <- "tur_t0m"
        names(turt0v) <- "tur_t0v"
        names(turMod) <- "turMod"
        names(turNLLslr) <- "tmslr_NLL"
        names(turAICslr) <- "tmslr_AIC"
        names(turBICslr) <- "tmslr_BIC"
        
        line_result <- t(c(temp_par, turMod, turNLLslr, turAICslr, turBICslr, 
                           turt0m, turt0v, temp_NLLsl, temp_AICsl, temp_BICsl,
                           temp_t0m,temp_t0v))
        total_result <- rbind(total_result, line_result)
    }
    total_result <- total_result[order(total_result[,5],
                                       total_result[,1],
                                       total_result[,2],
                                       total_result[,3],
                                       strtoi(total_result[,4])),]
    write.csv(total_result, file = filenm)
    total_result <- as.data.frame(total_result)
    list <- sapply(total_result, is.factor)
    list[length(par_pos)+1]<- F
    total_result[list] <- lapply(total_result[list], function(x) as.numeric(as.character(x)))
    return(total_result)
}

# Parameter setting -------------------------------------------------------

locf_set <- c("~1", "~t0", "~t", "~t+t0")
scaf_set <- c("~1")
par_pos <- data.frame(a1=1, a2=2, tt=3, xi=4, n=5)
t0_pos <- data.frame(t0tp1=12, t1tp1=14)
NLL_pos <- data.frame(t0tp0=7, t0tp1=8, t1tp0=9, t1tp1=10)
n_pos <- 5
k <- c(3, 5, 4, 6)
tuedep <- data.frame(t0tp0 = c(T, T),
                     t0tp1 = c(T, F),
                     t1tp0 = c(F, T),
                     t1tp1 = c(F, F))

# Simulation --------------------------------------------------------------

N = 1200 # excess 200 iteration for spare
n_set = c(30, 40, 50, 70, 90, 110, 150)
a1_set = c(-4, -2, -1, 1, 2, 4)
a2_set = c(-4, -2, -1, 1, 2, 4)
tt_set = c(0.25, 0.75)
xi_set = c(0.2, 0.1, 0, -0.1, -0.2)

simulation(a1_set, a2_set, tt_set, xi_set, n_set, N, locf_set, scaf_set)
simulation(0, a2_set, 0.75, xi_set, n_set, N, locf_set, scaf_set)
simulation(0, a2_set, 0.25, xi_set, n_set, N, locf_set, scaf_set)

test_0p50 <- report(..., par_pos, t0_pos,
                    NLL_pos, n_pos, k, tuedep)
test_0p25 <- report(..., par_pos, t0_pos,
                    NLL_pos, n_pos, k, tuedep)
test_0p75 <- report(..., par_pos, t0_pos,
                    NLL_pos, n_pos, k, tuedep)

# Plot --------------------------------------------------------------------

bpePlot <- function(total_result) {
    temp_plot1 = ggplot(total_result[total_result['turMod'] == "t0tp1",],
                        mapping = aes(x = as.factor(n), y = xi, fill= tur_t0v)) +
        facet_grid(. ~ a2, labeller = label_bquote(cols = beta[2]==.(x))) +
        geom_raster() + coord_cartesian(expand = FALSE) +
        scale_x_discrete(breaks=c(30,50,90,150), labels=c(30,50,90,150)) +
        scale_fill_gradientn(colours=c("#FFFFFF","#888888","#000000")) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        labs(x = TeX("$n$"), y = TeX("$\\xi$"), fill = TeX("$Var(\\tilde{t_b})$"),
             title = TeX("Breakpoint Estimation Variance of $\\mathrm{GEV}(0, 1), t_b = 0.5$"))
    ggsave(temp_plot1, file="tb_var_gev01.eps", device = "eps",
           width = 8, height = 2, dpi = 150, units = "in")
    
    temp_plot2 = ggplot(total_result[total_result['turMod'] == "t1tp1",],
                        mapping = aes(x = as.factor(n), y = xi, fill= tur_t0v)) +
        facet_grid(a1 ~ a2, labeller = label_bquote(cols = beta[2]==.(x),
                                                    rows = beta[1]==.(x))) +
        geom_raster() + coord_cartesian(expand = FALSE) +
        scale_x_discrete(breaks=c(30,50,90,150), labels=c(30,50,90,150)) +
        scale_fill_gradientn(colours=c("#FFFFFF","#888888","#000000")) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        labs(x = TeX("$n$"), y = TeX("$\\xi$"), fill = TeX("$Var(\\tilde{t_b})$"),
             title = TeX("Breakpoint Estimation Variance of $\\mathrm{GEV}(1, 1)$"))
    ggsave(temp_plot2, file="tb_var_gev11.eps", device = "eps",
           width = 8, height = 6, dpi = 150, units = "in")
    
    temp_plot3 = ggplot(total_result[total_result['turMod'] == "t0tp1",],
                        mapping = aes(x = as.factor(n), y = xi, fill= tur_t0m)) +
        facet_grid(. ~ a2, labeller = label_bquote(cols = beta[2]==.(x))) +
        geom_raster() + coord_cartesian(expand = FALSE) +
        scale_x_discrete(breaks=c(30,50,90,150), labels=c(30,50,90,150)) +
        scale_fill_gradientn(colours=c("#0000FF","#FFFFFF","#FF0000"),
                             breaks=c(0.45,0.5, 0.55), limits=c(0.45,0.55)) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        labs(x = TeX("$n$"), y = TeX("$\\xi$"), fill = TeX("$E(\\tilde{t_b})$"),
             title = TeX("Breakpoint Estimation Mean of $\\mathrm{GEV}(0, 1), t_b = 0.5$"))
    ggsave(temp_plot3, file="tb_mean_gev01.eps", device = "eps",
           width = 8, height = 2, dpi = 150, units = "in")
    
    temp_plot4 = ggplot(total_result[total_result['turMod'] == "t1tp1",],
                        mapping = aes(x = as.factor(n), y = xi, fill= tur_t0m)) +
        facet_grid(a1 ~ a2, labeller = label_bquote(cols = beta[2]==.(x),
                                                    rows = beta[1]==.(x))) +
        geom_raster() + coord_cartesian(expand = FALSE) +
        scale_x_discrete(breaks=c(30,50,90,150), labels=c(30,50,90,150)) +
        scale_fill_gradientn(colours=c("#0000FF","#FFFFFF","#FF0000"),
                             breaks=c(0.45,0.5, 0.55), limits=c(0.45,0.55)) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        labs(x = TeX("$n$"), y = TeX("$\\xi$"), fill = TeX("$E(\\tilde{t_b})$"),
             title = TeX("Breakpoint Estimation Mean of $\\mathrm{GEV}(1, 1)$"))
    ggsave(temp_plot4, file="tb_mean_gev11.eps", device = "eps",
           width = 8, height = 6, dpi = 150, units = "in")
}

bpePlotNC <- function(total_result, NC, c1, cc, c2) {
    temp_plot = ggplot(total_result[total_result['turMod'] == "t0tp1",],
                       mapping = aes(x = as.factor(n), y = xi, fill= tur_t0v)) +
        facet_grid(. ~ a2, labeller = label_bquote(cols = beta[2]==.(x))) +
        geom_raster() + coord_cartesian(expand = FALSE) +
        scale_x_discrete(breaks=c(30,50,90,150), labels=c(30,50,90,150)) +
        scale_fill_gradientn(colours=c("#FFFFFF","#888888","#000000")) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        labs(x = TeX("$n$"), y = TeX("$\\xi$"), fill = TeX("$Var(\\tilde{t_b})$"),
             title = TeX(paste0("Breakpoint Estimation Variance of $\\mathrm{GEV}(0, 1), t_b = ", NC,"$")))
    ggsave(temp_plot, file="tb_var_gev01.eps", device = "eps",
           width = 8, height = 2, dpi = 150, units = "in")
    
    temp_plot = ggplot(total_result[total_result['turMod'] == "t0tp1",],
                       mapping = aes(x = as.factor(n), y = xi, fill= tur_t0m)) +
        facet_grid(. ~ a2, labeller = label_bquote(cols = beta[2]==.(x))) +
        geom_raster() + coord_cartesian(expand = FALSE) +
        scale_x_discrete(breaks=c(30,50,90,150), labels=c(30,50,90,150)) +
        scale_fill_gradientn(colours=c("#0000FF","#FFFFFF","#FF0000"),
                             breaks=c(c1, cc, c2), limits=c(c1,c2)) +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        labs(x = TeX("$n$"), y = TeX("$\\xi$"), fill = TeX("$E(\\tilde{t_b})$"),
             title = TeX(paste0("Breakpoint Estimation Mean of $\\mathrm{GEV}(0, 1), t_b = ", NC,"$")))
    ggsave(temp_plot, file="tb_mean_gev01.eps", device = "eps",
           width = 8, height = 2, dpi = 150, units = "in")
    
}

pSelRate <- function(total_result){
    for(i in c(0.1,0.2,0,-0.1,-0.2)){
        title_t = paste0("Proper Selection Ratio of $\\mathrm{GEV}(1, 1),\\, \\xi = ",
                         i, "$", collapse = NULL)
        temp_plot = ggplot(total_result[total_result['turMod'] == "t1tp1" &
                                            total_result['xi'] == i,],
                           mapping = aes(x = n)) + 
            geom_line(aes(y = tmslr_NLL, colour = "NLL")) +
            geom_point(aes(y = tmslr_NLL, colour = "NLL", shape = "NLL")) +
            geom_line(aes(y = tmslr_AIC, colour = "AIC")) +
            geom_point(aes(y = tmslr_AIC, colour = "AIC", shape = "AIC")) +
            geom_line(aes(y = tmslr_BIC, colour = "BIC")) +
            geom_point(aes(y = tmslr_BIC, colour = "BIC", shape = "BIC")) +
            facet_grid(a1 ~ a2, labeller = label_bquote(cols = beta[2]==.(x),
                                                        rows = beta[1]==.(x))) +
            scale_x_continuous(breaks=c(30,70,110,150), labels=c(30,70,110,150)) +
            scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
            scale_color_discrete("Method") +
            scale_shape_manual("Method", values = c(16,17,15)) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = TeX("$n$"), y = "Proper model selection ratio", colour = "Method",
                 title = TeX(title_t))
        ggsave(temp_plot, file=paste0("Psr_gev11_xi", i,".eps"), device = "eps",
               width = 8, height = 6, dpi = 150, units = "in")
    }
    
    temp_plot = ggplot(total_result[total_result['turMod'] == "t1tp0",],
                       mapping = aes(x = n)) + 
        geom_line(aes(y = tmslr_NLL, colour = "NLL")) +
        geom_point(aes(y = tmslr_NLL, colour = "NLL", shape = "NLL")) +
        geom_line(aes(y = tmslr_AIC, colour = "AIC")) +
        geom_point(aes(y = tmslr_AIC, colour = "AIC", shape = "AIC")) +
        geom_line(aes(y = tmslr_BIC, colour = "BIC")) +
        geom_point(aes(y = tmslr_BIC, colour = "BIC", shape = "BIC")) +
        facet_grid(xi ~ a1, labeller = label_bquote(cols = beta[1]==.(x),
                                                    rows = xi==.(x))) +
        scale_x_continuous(breaks=c(30,70,110,150), labels=c(30,70,110,150)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        scale_color_discrete("Method") +
        scale_shape_manual("Method", values = c(16,17,15)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = TeX("$n$"), y = "Proper model selection ratio", colour = "Method",
             title = TeX("Proper Selection Ratio of $\\mathrm{GEV}(1, 0)$"))
    ggsave(temp_plot, file="Psr_gev10.eps", device = "eps",
           width = 8, height = 6, dpi = 150, units = "in")
    
    temp_plot = ggplot(total_result[total_result['turMod'] == "t0tp1",],
                       mapping = aes(x = n)) + 
        geom_line(aes(y = tmslr_NLL, colour = "NLL")) +
        geom_point(aes(y = tmslr_NLL, colour = "NLL", shape = "NLL")) +
        geom_line(aes(y = tmslr_AIC, colour = "AIC")) +
        geom_point(aes(y = tmslr_AIC, colour = "AIC", shape = "AIC")) +
        geom_line(aes(y = tmslr_BIC, colour = "BIC")) +
        geom_point(aes(y = tmslr_BIC, colour = "BIC", shape = "BIC")) +
        facet_grid(xi ~ a2, labeller = label_bquote(cols = beta[2]==.(x),
                                                    rows = xi==.(x))) +
        scale_x_continuous(breaks=c(30,70,110,150), labels=c(30,70,110,150)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        scale_color_discrete("Method") +
        scale_shape_manual("Method", values = c(16,17,15)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = TeX("$n$"), y = "Proper model selection ratio", colour = "Method",
             title = TeX("Proper Selection Ratio of $\\mathrm{GEV}(0, 1)$"))
    ggsave(temp_plot, file="Psr_gev01.eps", device = "eps",
           width = 8, height = 6, dpi = 150, units = "in")
    
    temp_plot = ggplot(total_result[total_result['turMod'] == "t0tp0",],
                       mapping = aes(x = n)) + 
        geom_line(aes(y = tmslr_NLL, colour = "NLL")) +
        geom_point(aes(y = tmslr_NLL, colour = "NLL", shape = "NLL")) +
        geom_line(aes(y = tmslr_AIC, colour = "AIC")) +
        geom_point(aes(y = tmslr_AIC, colour = "AIC", shape = "AIC")) +
        geom_line(aes(y = tmslr_BIC, colour = "BIC")) +
        geom_point(aes(y = tmslr_BIC, colour = "BIC", shape = "BIC")) +
        facet_grid(. ~ xi, labeller = label_bquote(cols = xi==.(x))) +
        scale_x_continuous(breaks=c(30,70,110,150), labels=c(30,70,110,150)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        scale_color_discrete("Method") +
        scale_shape_manual("Method", values = c(16,17,15)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = TeX("$n$"), y = "Proper model selection ratio", colour = "Method",
             title = TeX("Proper Selection Ratio of $\\mathrm{GEV}(0, 0)$"))
    ggsave(temp_plot, file="Psr_gev00.eps", device = "eps",
           width = 8, height = 2, dpi = 150, units = "in")
    
    a1_set = c(-4, -2, -1, 0, 1, 2, 4)
    
    for(i in a1_set[a1_set!=0]){
        title_t = paste0("Proper Selection Ratio of $\\mathrm{GEV}(1, 1),\\, \\beta_1 = ",
                         i, "$", collapse = NULL)
        temp_plot = ggplot(total_result[total_result['turMod'] == "t1tp1" &
                                            total_result['a1'] == i,],
                           mapping = aes(x = n)) + 
            geom_line(aes(y = tmslr_NLL, colour = "NLL")) +
            geom_point(aes(y = tmslr_NLL, colour = "NLL", shape = "NLL")) +
            geom_line(aes(y = tmslr_AIC, colour = "AIC")) +
            geom_point(aes(y = tmslr_AIC, colour = "AIC", shape = "AIC")) +
            geom_line(aes(y = tmslr_BIC, colour = "BIC")) +
            geom_point(aes(y = tmslr_BIC, colour = "BIC", shape = "BIC")) +
            facet_grid(xi ~ a2, labeller = label_bquote(cols = beta[2]==.(x),
                                                        rows = xi ==.(x))) +
            scale_x_continuous(breaks=c(30,70,110,150), labels=c(30,70,110,150)) +
            scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
            scale_color_discrete("Method") +
            scale_shape_manual("Method", values = c(16,17,15)) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = TeX("$n$"), y = "Proper model selection ratio", colour = "Method",
                 title = TeX(title_t))
        ggsave(temp_plot, file=paste0("Psr_gev11_beta1_", i,".eps"), device = "eps",
               width = 8, height = 6, dpi = 150, units = "in")
    }
}

pSelRateNC <- function(total_result, NC){
    temp_plot = ggplot(total_result[total_result['turMod'] == "t0tp1",],
                       mapping = aes(x = n)) + 
        geom_line(aes(y = tmslr_NLL, colour = "NLL")) +
        geom_point(aes(y = tmslr_NLL, colour = "NLL", shape = "NLL")) +
        geom_line(aes(y = tmslr_AIC, colour = "AIC")) +
        geom_point(aes(y = tmslr_AIC, colour = "AIC", shape = "AIC")) +
        geom_line(aes(y = tmslr_BIC, colour = "BIC")) +
        geom_point(aes(y = tmslr_BIC, colour = "BIC", shape = "BIC")) +
        facet_grid(xi ~ a2, labeller = label_bquote(cols = beta[2]==.(x),
                                                    rows = xi==.(x))) +
        scale_x_continuous(breaks=c(30,70,110,150), labels=c(30,70,110,150)) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
        scale_color_discrete("Method") +
        scale_shape_manual("Method", values = c(16,17,15)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = TeX("$n$"), y = "Proper model selection ratio", colour = "Method",
             title = TeX(paste0("Proper Selection Ratio of $\\mathrm{GEV}(0, 1), t_b = ", NC,"$")))
    ggsave(temp_plot, file="Psr_gev01.eps", device = "eps",
           width = 8, height = 6, dpi = 150, units = "in")
}

mSelRate <- function(total_result){
    tr_temp <- gather(total_result,
                      'NLL_t0tp0', 'NLL_t0tp1', 'NLL_t1tp0', 'NLL_t1tp1',
                      'AIC_t0tp0', 'AIC_t0tp1', 'AIC_t1tp0', 'AIC_t1tp1',
                      'BIC_t0tp0', 'BIC_t0tp1', 'BIC_t1tp0', 'BIC_t1tp1',
                      key = "Method_Model", value = "Sl_num")
    tr_temp <- separate(tr_temp, Method_Model, sep = "_",
                        into = c("Sl_Meth", "Sl_Model"))
    tr_temp$Sl_Model[tr_temp["Sl_Model"] == "t0tp0"] = "GEV(0,0)"
    tr_temp$Sl_Model[tr_temp["Sl_Model"] == "t1tp0"] = "GEV(1,0)"
    tr_temp$Sl_Model[tr_temp["Sl_Model"] == "t0tp1"] = "GEV(0,1)"
    tr_temp$Sl_Model[tr_temp["Sl_Model"] == "t1tp1"] = "GEV(1,1)"
    tr_temp$Sl_Model <- as.factor(tr_temp$Sl_Model)
    tr_temp$Sl_Meth <- as.factor(tr_temp$Sl_Meth)
    tr_temp$n <- as.factor(tr_temp$n)
    tr_temp$Sl_num <- tr_temp$Sl_num/1000
    
    n_set = c(30, 40, 50, 70, 90, 110, 150)
    a1_set = c(-4, -2, -1, 0, 1, 2, 4)
    a2_set = c(-4, -2, -1, 0, 1, 2, 4)
    xi_set = c(-0.2, -0.1, 0, 0.1, 0.2)
    
    for(a1_i in a1_set){
        for(a2_i in a2_set){
            for(xi_i in xi_set){
                if(a1_i != 0){modt_t <- 1}else{modt_t <- 0}
                if(a2_i != 0){modtp_t <- 1}else{modtp_t <- 0}
                "Model Selection Ratio of $\\mathrm{GEV}(0, 0)$"
                title_t = paste0("Model selection ratio of $\\mathrm{GEV}(",
                                 modt_t, ",", modtp_t, "),\\, \\xi = ",
                                 xi_i, ",\\, \\beta_1 = ", a1_i,
                                 ",\\, \\beta_2 = ", a2_i,"$", collapse = NULL)
                fn_t = paste0("Msr_gev", modt_t, modtp_t, "_xi", xi_i, "_a1", a1_i,
                              "_a2", a2_i, ".eps", collapse = NULL)
                
                temp_plot = ggplot(tr_temp[tr_temp['a1'] == a1_i&
                                               tr_temp['a2'] == a2_i&
                                               tr_temp['xi'] == xi_i,],
                                   mapping = aes(x = n, y = Sl_num, fill = Sl_Model)) + 
                    geom_bar(position="dodge", stat="identity") + 
                    facet_grid(.~ Sl_Meth) + 
                    scale_x_discrete(breaks=n_set, labels=n_set) +
                    scale_y_continuous(labels = scales::percent, limits = c(0,1)) + 
                    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs(x = TeX("$n$"), y = "Model selection ratio", fill = "Model",
                         title = TeX(title_t))
                
                ggsave(temp_plot, file=fn_t, device = "eps",
                       width = 8, height = 3, dpi = 150, units = "in")
            }
        }
    }
}


bpePlot(test_0p50)
pSelRate(test_0p50)
mSelRate(test_0p50)

bpePlotNC(test_0p25, "0.25", 0,0.25,0.5)
bpePlotNC(test_0p75, "0.75", 0.5,0.75,1)

pSelRateNC(test_0p25, "0.25")
pSelRateNC(test_0p75, "0.75")

# 以下 GEV01_tbMLE.csv, GEV11_tbMLE.csv 為整理單模型模擬結果之資料表
GEV01_tbMLE <- read_csv("GEV01_tbMLE.csv", 
    col_types = cols(beta_1 = col_factor(levels = c("0")), 
        beta_2 = col_factor(levels = c("-4", "-2", "-1", "1", "2", "4")), 
        n = col_factor(levels = c("30", "40", "50", "70", "90", "110", "150")), 
        tb = col_factor(levels = c("0.25", "0.5", "0.75")),
        xi = col_factor(levels = c("-0.2", "-0.1", "0", "0.1", "0.2"))))

temp_1 <- ggplot(data = GEV01_tbMLE) +
    aes(x = tb, y = tb_mean, fill = tb) +
    geom_violin(scale = "count", adjust = 1) +
    labs(title = TeX("Mean of Estimate Breakpoint Position in $\\mathrm{GEV}(0,1)$"),
         x = "Setting position of breakpoint",
         y = "Mean of estimate position") +
    theme_minimal() +
    theme(legend.position = 'none')

ggsave(temp_plot, file="MEBP_GEV01.eps", device = "eps",
       width = 8, height = 4, dpi = 150, units = "in")

temp_2 <- ggplot(data = GEV01_tbMLE) +
    aes(x = tb, y = tb_var, fill = tb) +
    geom_violin(scale = "count", adjust = 1) +
    labs(title = TeX("Variance of Estimate Breakpoint Position in $\\mathrm{GEV}(0,1)$"),
         x = "Setting position of breakpoint",
         y = "Variance of estimate position") +
    theme_minimal() +
    theme(legend.position = 'none')

temp_plot <-plots(temp_1, temp_2)

ggsave(temp_plot, file="MEBP_GEV01.eps", device = "eps",
       width = 8, height = 8, dpi = 150, units = "in")

GEV11_tbMLE <- read_csv("GEV11_tbMLE.csv", 
    col_types = cols(beta_1 = col_factor(levels = c("-4", "-2", "-1", "1", "2", "4")),
        beta_2 = col_factor(levels = c("-4", "-2", "-1", "1", "2", "4")),
        n = col_factor(levels = c("30", "40", "50", "70", "90", "110", "150")), 
        tb = col_factor(levels = c("0.5")), 
        xi = col_factor(levels = c("-0.2", "-0.1", "0", "0.1", "0.2"))))

temp_3 <- ggplot(data = GEV11_tbMLE) +
    aes(x = beta_2, y = tb_mean, fill = beta_2) +
    geom_violin(scale = "area", adjust = 1) +
    labs(title = TeX("Mean of Estimate Breakpoint Position in $\\mathrm{GEV}(1,1)$"),
         x = "Setting gap size of breakpoint",
         y = "Mean of estimate position") +
    theme_minimal() +
    theme(legend.position = 'none')

temp_4 <- ggplot(data = GEV11_tbMLE) +
    aes(x = beta_2, y = tb_var, fill = beta_2) +
    geom_violin(scale = "area", adjust = 1) +
    labs(title = TeX("Variance of Estimate Breakpoint Position in $\\mathrm{GEV}(1,1)$"),
         x = "Setting gap size of breakpoint",
         y = "Mean of estimate position") +
    theme_minimal() +
    theme(legend.position = 'none')

temp_plot_2 <-plots(temp_3, temp_4)

ggsave(temp_plot_2, file="MEBP_GEV11.eps", device = "eps",
       width = 8, height = 8, dpi = 150, units = "in")
