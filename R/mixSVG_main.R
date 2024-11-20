mixSVG_main <- function (y, X, s_trans, pat_idx, pat_name, perm_sample, libsize, vtest_zero_prop, method='davies'){
    vtest = (mean(y == 0) < vtest_zero_prop)
    model_init = glm(y ~ X - 1 + offset(log(libsize)), family = poisson)
    model0 = fit_glmm(y, X, model_init, libsize)
    par = model0$par[1, ]
    beta = par[1:ncol(X)]
    tau = par['tau']
    w = model0$w
    vw = model0$vw
    mu =  model0$mu
    XVivX_iv = solve(t(X/vw)%*%X)
    res = (w - X %*% beta)/vw
    res2 = res^2
    res2_perm = matrix(res2[perm_sample], nrow = nrow(perm_sample))


    test_func = function(i_pat) {
        s1 = s_trans[, (2 * i_pat - 1)]
        s2 = s_trans[, (2 * i_pat)]
        s = cbind(s1, s2)
        s1_sq = s1^2
        s2_sq = s2^2
        Tb = c(sum(res * s1), sum(res * s2))
        ZivX = t(s/vw) %*% X
        Vb = t(s/vw) %*% s - ZivX %*% XVivX_iv %*% t(ZivX)
        Tb = t(Tb) %*% solve(Vb) %*% t(t(Tb))
        pval_b = pchisq(Tb, 2, lower.tail = F)
      
        if (vtest) {
            s_sq = s1_sq + s2_sq
            Tv = sum(res2 * s_sq)
      
            PK = diag(s_sq/vw) - (X/vw)%*%XVivX_iv%*%t(X/vw*s_sq)
            lambdas = eigen(PK)$values
            lambdas = lambdas[!is.complex(lambdas)]
            if(method=="davies"){
              pval_v = davies(Tv,lambdas)$Qq
              if(pval_v == 0 | pval_v < 0){
                pval_v = liumod(Tv,lambdas)
              }
            }else if(method=="liu"){
              pval_v = liu(Tv,lambdas) 
            }else if(method=="liumod"){
              pval_v = liumod(Tv,lambdas)
            }
            
            
            pval = c(pval_b, pval_v)
            pval[which(pval == 0)] <- 5.55e-17
            pval[which((1 - pval) < 0.001)] <- 0.99
            T_omn = mean(tan(pi * (0.5 - pval)))
            pval = 1 - pcauchy(T_omn)
        }
        else {
            pval_v = 1
            pval = pval_b
            Tv =  k = df = 0
        }
        return(c(pval, pval_b, pval_v))
    }
    pval_pat = t(apply(t(pat_idx), 2, FUN = test_func))
    colnames(pval_pat) = c("pval_omn", "pval_b", "pval_v")
    pval = pval_pat[, "pval_omn"]
    pval[which(pval == 0)] <- 5.55e-17
    pval[which((1 - pval) < 0.001)] <- 0.99
    T_final = mean(tan(pi * (0.5 - pval)))
    pval = 1 - pcauchy(T_final)
    out = list(model0 = par, pval = pval, pval_pat = pval_pat, 
        pattern = pat_name, res = res, vw = vw, mu = mu, w = w  )
    return(out)
}
