library(bayesianVARs)
setwd("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code")
load("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/StartingPoint.RData")
rm(list = setdiff(ls(), c("macro", "data", "cc.var.mn", "cc.var.mn.sv")))


all_models <- function(p){
  cc.var.chol.Nphi.Nchol <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                               prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "normal"))
  cc.var.chol.Nphi.HSchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "hs"))
  cc.var.chol.Nphi.NGchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "ng"))
  cc.var.chol.Nphi.SSVSchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "ssvs"))


  cc.var.chol.NGphi.Nchol <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "normal"))
  cc.var.chol.NGphi.HSchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "hs"))
  cc.var.chol.NGphi.NGchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "ng"))
  cc.var.chol.NGphi.SSVSchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "ssvs"))

  cc.var.chol.SSVSphi.Nchol <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                  prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "normal"))
  cc.var.chol.SSVSphi.HSchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "hs"))
  cc.var.chol.SSVSphi.NGchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "ng"))
  cc.var.chol.SSVSphi.SSVSchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "ssvs"))

  cc.var.chol.HSphi.Nchol <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "normal"))
  cc.var.chol.HSphi.HSchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "hs"))
  cc.var.chol.HSphi.NGchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "ng"))
  cc.var.chol.HSphi.SSVSchol <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="cholesky", cholesky_U_prior = "ssvs"))
  #7 Factor w/out minimum on factor


  cc.var.fsv.Nphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                              prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                              prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.Nphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                               prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.Nphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                               prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.Nphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))


  cc.var.fsv.NGphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                               prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.NGphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.NGphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.NGphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))

  cc.var.fsv.SSVSphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                 prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.SSVSphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.SSVSphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.SSVSphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                    prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))

  cc.var.fsv.HSphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                               prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.HSphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.HSphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.HSphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))

  cc.var.fsv.Nphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                              prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                              prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.Nphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                               prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.Nphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                               prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.Nphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))


  cc.var.fsv.NGphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                               prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.NGphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.NGphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.NGphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))

  cc.var.fsv.SSVSphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                 prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.SSVSphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.SSVSphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.SSVSphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                    prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))

  cc.var.fsv.HSphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                               prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.HSphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.HSphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  cc.var.fsv.HSphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=7))
  #17 factor w/ out minimum on factor

  cc.var.fsv.Nphi.Nfsv17 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                               prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                               prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.Nphi.HSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.Nphi.NGfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.Nphi.SSVSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))


  cc.var.fsv.NGphi.Nfsv17 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.NGphi.HSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.NGphi.NGfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.NGphi.SSVSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))

  cc.var.fsv.SSVSphi.Nfsv17 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                  prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.SSVSphi.HSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.SSVSphi.NGfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.SSVSphi.SSVSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))

  cc.var.fsv.HSphi.Nfsv17 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.HSphi.HSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.HSphi.NGfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                 prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))
  cc.var.fsv.HSphi.SSVSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_factors=17))


  #7 Factor w/ factor minimum

  cc.var.minFSV.Nphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                 prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.Nphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.Nphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.Nphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                    prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))


  cc.var.minFSV.NGphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                  prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.NGphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.NGphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.NGphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))

  cc.var.minFSV.SSVSphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                    prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.SSVSphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.SSVSphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.SSVSphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                       prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                       prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))

  cc.var.minFSV.HSphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                  prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.HSphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.HSphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.HSphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))

  cc.var.minFSV.Nphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                 prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                 prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.Nphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.Nphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                  prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.Nphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                    prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))


  cc.var.minFSV.NGphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                  prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.NGphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.NGphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.NGphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))

  cc.var.minFSV.SSVSphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                    prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.SSVSphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.SSVSphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.SSVSphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                       prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                       prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))

  cc.var.minFSV.HSphi.Nfsv7 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                  prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.HSphi.HSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.HSphi.NGfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  cc.var.minFSV.HSphi.SSVSfsv7 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=7))
  #17 factor w/ factor minimum

  cc.var.minFSV.Nphi.Nfsv17 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                  prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                  prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.Nphi.HSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.Nphi.NGfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                   prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.Nphi.SSVSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                     prior_phi = specify_prior_phi(data=data, prior="normal", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))


  cc.var.minFSV.NGphi.Nfsv17 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                   prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.NGphi.HSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                    prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.NGphi.NGfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                    prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.NGphi.SSVSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                      prior_phi = specify_prior_phi(data=data, prior="ng", lags=p),
                                                      prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))

  cc.var.minFSV.SSVSphi.Nfsv17 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                     prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                     prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.SSVSphi.HSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                      prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                      prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.SSVSphi.NGfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                      prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                      prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.SSVSphi.SSVSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                        prior_phi = specify_prior_phi(data=data, prior="SSVS", lags=p),
                                                        prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))

  cc.var.minFSV.HSphi.Nfsv17 <- bayesianVARs::bvar(data, p, sv_keep = "all",
                                                   prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                   prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.HSphi.HSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                    prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.HSphi.NGfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                    prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                    prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))
  cc.var.minFSV.HSphi.SSVSfsv17 <- bayesianVARs::bvar(data, p, sv_keep="all",
                                                      prior_phi = specify_prior_phi(data=data, prior="HS", lags=p),
                                                      prior_sigma = specify_prior_sigma(data=data, type="factor", factor_restrict = "none", factor_facloadtol=1e-18, factor_factors=17))

  save(list = ls(), file = paste0("CC_VARs_bayesianVARs_fsv_chol_lag", p, ".RData"), envir = environment())
}

all_models(4)








