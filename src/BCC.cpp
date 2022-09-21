#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
#include <math.h>

// #define DEBUG
//#define NORAND

// - AlignCllusters function - adapted from the codes of BCC original paper (Lock and Dunson 2013)
arma::rowvec AlignClusters (arma::rowvec Z1, arma::rowvec Z2, std::string type = "vec") {
  // if (type.compare("vec") == 0) {
    arma::rowvec uniq1 = unique(Z1);
    arma::rowvec uniq2 = unique(Z2);
    for (int k = 0; k < uniq1.size(); k++) {
      double Max = arma::sum(Z1==k && Z2==k)/(.01+arma::sum(Z2==k)+arma::sum(Z1==k));
      for(int tempk = 0; tempk < uniq2.size(); tempk++) {
        if(arma::sum(Z1==k && Z2==tempk)/(.01+arma::sum(Z2==tempk)+arma::sum(Z1==k)) > Max) {
          Max                 = arma::sum(Z1==k && Z2==tempk)/(.01+arma::sum(Z2==tempk)+arma::sum(Z1==k));
          arma::urowvec dummy = Z2==k;
          Z2.elem(Z2==tempk) *= 0;
          Z2.elem(Z2==tempk) += k;
          Z2.elem(dummy)     *= 0;
          Z2.elem(dummy)     += tempk;
        }
      }
    }
  // currently not being used
  // } else if(type.compare("mat") == 0) {
  //   for(int k = 0; k < dim(Z1)[2]; k++) {
  //     for(tempk in  1:dim(Z2)[2]) {
  //       Max             = sum(Z1==Z2)
  //       Z2dummy         = Z2
  //       Z2dummy[,k]     = Z2[,tempk]
  //       Z2dummy[,tempk] = Z2[,k]
  //       if(sum(Z1==Z2dummy)>Max)
  //         Z2 = Z2dummy
  //     }
  //   }
  // }
  return Z2;
}

// [[Rcpp::export()]]
Rcpp::List BCC (
    Rcpp::List dat, int R,
    Rcpp::List id, arma::umat n_obs, int N,
    int num_cluster,
    std::vector<std::string> dist,
    bool alpha_common,
    bool sigma_sq_e_common,
    bool align_clusters,
    arma::vec p,
    arma::vec q,
    bool sig_var,

    arma::vec ppi,
    arma::vec alpha,
    arma::vec zz,
    arma::mat zz_local,
    std::vector<arma::mat> gamma,
    std::vector<arma::vec> sigma_sq_e,
    std::vector<arma::vec> phi,
    std::vector<arma::cube> sigma_sq_u,
    std::vector<arma::mat> beta,

    double delta,
    double a_star,
    double b_star,
    double aa0,
    double bb0,
    arma::mat a0,
    arma::mat b0,
    std::vector<arma::mat> v0,
    std::vector<arma::cube> V0,
    double cc0,
    double dd0,
    arma::mat c0,
    arma::mat d0,
    double rr0,
    double RR0,
    double ww0,
    double vv0,
    arma::mat lambda0,
    std::vector<arma::cube> Lambda0,
    
             Rcpp::RObject LOG_LIK_ITER,
             Rcpp::RObject PPI,
             Rcpp::RObject ZZ,
             Rcpp::RObject ALPHA,
    std::vector<arma::mat> ZZ_LOCAL,
    std::vector<Rcpp::RObject> GA,
    std::vector<arma::mat> GA_ACCEPT,
    std::vector<arma::cube> THETA,
    std::vector<arma::mat> THETA_ACCEPT,
    std::vector<Rcpp::RObject> SIGMA_SQ_U,
    std::vector<arma::mat> SIGMA_SQ_E,
    std::vector<arma::cube> T_LOCAL,
             Rcpp::RObject T,

    bool adaptive_tuning,
    double tuning_freq,
    arma::mat c_gamma_tuning,
    arma::vec c_beta_tuning,

    int burn_in,
    int thin,
    int per,
    int max_iter,
    int seed_initial
  ) {
  Rcpp::Rcout << "Cpp hi" << std::endl ;
#ifdef NORAND
  Rcpp::Function set_seed("set.seed");
#endif
  Rcpp::Environment mclust = Rcpp::Environment::namespace_env("mclust");
  Rcpp::Function dmvnorm = mclust["dmvnorm"];
  Rcpp::Function dpois("dpois");
  Rcpp::Function dbinom("dbinom");

  Rcpp::Function rbind("rbind");
  Rcpp::Environment abindenv = Rcpp::Environment::namespace_env("abind");
  Rcpp::Function abind = abindenv["abind"];

  Rcpp::Function apply("apply");
  Rcpp::Function mean("mean");

  int count = 0;
  int iter = 0;
  Rcpp::List rst = Rcpp::List::create();

  while(1) {
    //--------------------------------------------------------------#
    // Sample outcome-specific cluster membership
    //--------------------------------------------------------------#
#ifdef DEBUG
    Rcpp::Rcout << "Sample outcome-specific cluster membership" << std::endl;
#endif
    arma::mat ppt;
    arma::cube pt;
    if (num_cluster > 1) {
      pt = arma::cube(N, num_cluster, R);
      arma::cube  v(N, num_cluster, R);
      arma::cube  f(N, num_cluster, R);
      for (int r = 0; r < R; r++) {
        arma::vec       y = ((Rcpp::DataFrame)dat(r))["y"];
        arma::vec       t = ((Rcpp::DataFrame)dat(r))["time"];
        arma::vec      id = ((Rcpp::DataFrame)dat(r))["id"];
        arma::uvec n_obss =             n_obs.col(r);

        for (int i = 0; i < N; i++) {
          arma::vec  in_i = arma::conv_to<arma::vec>::from(id==i+1); // c++ index begins from 0
          arma::uvec in_i_idx = find(in_i);
          arma::vec  t_in_i = t.elem(in_i_idx);
          arma::vec one(t_in_i.size(), arma::fill::ones);
          arma::mat x_in_i = arma::join_rows(one, t_in_i, pow(t_in_i, 2), pow(t_in_i, 3));
          arma::mat y_in_i = y.elem(in_i_idx);
          arma::mat x  = x_in_i.cols(0, p[r]-1);
          arma::mat Z  = x_in_i.cols(0, q[r]-1);

          for (int k = 0; k < num_cluster; k++) {
            arma::rowvec gamma_(gamma[r].row(k).subvec(0, p[r]-1));
            arma::rowvec eta = gamma_ * x.t() + beta[r].row(i) * Z.t();
            arma::rowvec mu;

            if (dist[r].compare("gaussian") == 0) {
              mu = eta;

              arma::mat sigma(n_obss[i], n_obss[i], arma::fill::zeros);
              sigma.diag() += sigma_sq_e[r][k];

              f(i,k,r) = *REAL(dmvnorm(
                y_in_i.t(),
                Rcpp::_["mean"]  = mu,
                Rcpp::_["sigma"] = sigma
              ));
            } else if (dist[r].compare("poisson") == 0) {
              mu = exp(eta); // TODO is this right?
              f(i,k,r) = prod(Rcpp::as<arma::vec>(dpois(
                y_in_i.t(),
                Rcpp::_["lambda"] = mu
              )));
            } else if (dist[r].compare("binomial") == 0) {
              mu = exp(eta)/(1+exp(eta)); // TODO is this right?
              f(i,k,r) = prod(Rcpp::as<arma::vec>(dbinom(
                y_in_i.t(),
                Rcpp::_["size"] = 1,
                Rcpp::_["prob"] = mu // mean = np && n = 1 ==> p = mean
              )));
            }
            v(i,k,r)  = (zz[i]==(k+1)) ? alpha[r] : (1-alpha[r])/(num_cluster-1);
            pt(i,k,r) = v(i,k,r)*f(i,k,r);
            // pt(i,k,r) = std::max(pt(i,k,r),1e-6);
          }    // proportional to
          arma::ivec dist(num_cluster);
          arma::rowvec ptc = pt.slice(r).row(i);
          arma::rowvec ptc_normed = ptc / arma::sum(ptc); // unlike R's internal rmultinom, Rcpp's rmultinom does not internally normalize distribution
#ifdef NORAND
          set_seed(seed_initial + iter + r+1 + i+1);          
#endif
          rmultinom(1, ptc_normed.begin(), num_cluster, dist.begin());

          zz_local(r, i) = dist.index_max() + 1; // R index starts from 1
        }
        if (align_clusters) {
          zz_local.row(r) = AlignClusters(zz.t(), zz_local.row(r), "vec");
        }
      }
#ifdef DEBUG
      rst.push_back(zz_local ,"zz.local");
#endif

      // mean(zz.local[[3]]==simdata$cluster.local[[3]])
      // table(zz.local[[3]],simdata$cluster.local[[3]]) 

      //--------------------------------------------------------------#
      // Sample zz: overall clustering
      //--------------------------------------------------------------#
#ifdef DEBUG
      Rcpp::Rcout << "Sample zz: overall clustering" << std::endl;
#endif
      ppt = arma::mat(N, num_cluster);
      arma::vec v2(R);
      for (int i = 0; i < N; i++) {
        for (int k = 0; k < num_cluster; k++) {
          for (int r = 0; r < R; r++) {
            v2[r] = (zz_local(r,i)==k+1) ? alpha[r] : (1-alpha[r])/(num_cluster-1);
          }
          ppt(i,k) = ppi[k] * prod(v2);
        }
        Rcpp::IntegerVector dist(num_cluster);
        arma::rowvec pptc = ppt.row(i);
        arma::rowvec pptc_normed = pptc / arma::sum(pptc); // unlike R's internal rmultinom, Rcpp's rmultinom does not internally normalize distribution
#ifdef NORAND
        set_seed(seed_initial + iter + i+1);
#endif
        rmultinom(1, pptc_normed.begin(), num_cluster, dist.begin());

        Rcpp::IntegerVector s = Rcpp::seq(0, dist.size()-1);
        Rcpp::LogicalVector compare = dist==1;
        Rcpp::IntegerVector pos = s[compare];
        zz[i] = pos[0] + 1; // R index starts from 1
      }
    }
#ifdef DEBUG
    rst.push_back(zz ,"zz");
#endif

    //--------------------------------------------------------------#
    //  Sample alpha
    //--------------------------------------------------------------#
#ifdef DEBUG
    Rcpp::Rcout << "Sample alpha" << std::endl;
#endif
#ifdef NORAND
    set_seed(seed_initial+iter);
#endif
    if (R > 1 && num_cluster > 1) {
        Rcpp::Environment truncdist = Rcpp::Environment::namespace_env("truncdist");
        Rcpp::Function rtrunc = truncdist["rtrunc"];
        Rcpp::NumericVector tau = Rcpp::no_init(zz_local.n_rows);
        for (int i = 0; i < zz_local.n_rows; i++) {
          arma::rowvec row = zz_local.row(i);
          tau(i) = arma::sum(row.t() == zz);
        }
        if (alpha_common == 1) {
          int tau_sum = sum(tau);
          alpha   = arma::vec(R,
                              *REAL(rtrunc(1,
                                          Rcpp::_["spec"]   = "beta",
                                          Rcpp::_["a"]      = 1.0/num_cluster,
                                          Rcpp::_["b"]      = 1,
                                          Rcpp::_["shape1"] = a_star       + tau_sum,
                                          Rcpp::_["shape2"] = b_star + N*R - tau_sum
                              )));
        } else {
          for (int r = 0; r < R; r++) {
            alpha[r] = *REAL(rtrunc(1,
                                    Rcpp::_["spec"]   = "beta",
                                    Rcpp::_["a"]      = 1.0/num_cluster,
                                    Rcpp::_["b"]      = 1,
                                    Rcpp::_["shape1"] = a_star     + tau[r],
                                    Rcpp::_["shape2"] = b_star + N - tau[r]
                                  ));
          }
        }
    } else {
      alpha = 1;
    }
#ifdef DEBUG
    rst.push_back(alpha ,"alpha");
#endif
      
    //--------------------------------------------------------------#
    //  Sample pi
    //--------------------------------------------------------------#
#ifdef DEBUG
    Rcpp::Rcout << "Sample pi" << std::endl;
#endif
    if (num_cluster > 1) {
      Rcpp::Environment MCMCpack = Rcpp::Environment::namespace_env("MCMCpack");
      Rcpp::Function rdirichlet = MCMCpack["rdirichlet"];
      arma::vec rho(num_cluster);
      for (int k = 0; k < num_cluster; k++) { rho[k] = arma::sum(zz==(k+1)); }
#ifdef NORAND
      set_seed(seed_initial+iter);
#endif
      ppi = Rcpp::as<arma::vec>(rdirichlet(1, delta + rho));
    } else {
      ppi = 1;
    }
#ifdef DEBUG
    rst.push_back(ppi ,"ppi");
#endif

    //--------------------------------------------------------------#
    // Sample Fixed Effect Coefficients (GAMMA) Using MH Alogrithm
    //--------------------------------------------------------------#
#ifdef DEBUG
    Rcpp::Rcout << "Sample Fixed Effect Coefficients (GAMMA) Using MH Alogrithm" << std::endl;
#endif
    std::vector<arma::cube> H(R);
    std::vector<arma::mat>  G(R);
    for  (int r = 0; r < R; r++) {
      H[r] = arma::cube(p[r], p[r], num_cluster);
      G[r] = arma::mat(num_cluster, p[r]);
      arma::vec       y = ((Rcpp::DataFrame)dat(r))["y"];
      arma::vec       t = ((Rcpp::DataFrame)dat(r))["time"];
      arma::vec      id = ((Rcpp::DataFrame)dat(r))["id"];
      arma::uvec n_obss =             n_obs.col(r);

      for (int i = 0; i < N; i++) {
        arma::vec  in_i = arma::conv_to<arma::vec>::from(id==i+1); // c++ index begins from 0
        arma::uvec in_i_idx = find(in_i);
        arma::vec  t_in_i = t.elem(in_i_idx);
        arma::vec one(t_in_i.size(), arma::fill::ones);
        arma::mat x_in_i = arma::join_rows(one, t_in_i, pow(t_in_i, 2), pow(t_in_i, 3));
        arma::mat y_in_i = y.elem(in_i_idx);
        arma::mat x  = x_in_i.cols(0, p[r]-1);
        arma::mat Z  = x_in_i.cols(0, q[r]-1);

        for(int k = 0; k < num_cluster; k++) {
          arma::rowvec gamma_(gamma[r].row(k).subvec(0, p[r]-1));
          arma::rowvec eta = gamma_ * x.t() + beta[r].row(i) * Z.t();
          arma::rowvec mu;
          arma::mat kappa(x.n_rows, x.n_rows);
          
          if (dist[r].compare("gaussian") == 0) {
            mu = eta;
            kappa.diag() += sigma_sq_e[r][k];
          } else if (dist[r].compare("poisson") == 0) {
            mu = exp(eta); // TODO is this right?
            kappa.diag() += exp(eta);
          } else if (dist[r].compare("binomial") == 0) {
            mu = exp(eta)/(1+exp(eta)); // TODO is this right?
            arma::rowvec p = exp(eta)/(1+exp(eta));
            kappa.diag() += p % (1-p); // element-wise product
          }

          H[r].slice(k) +=  pow(phi[r][k],-2) * (zz_local(r,i)==(k+1)) * x.t() * kappa * x;
          G[r].row(k)   += (pow(phi[r][k],-1) * (zz_local(r,i)==(k+1)) * x.t() * (y_in_i.t() - mu).t()).t();
        }
      }
    }
#ifdef DEBUG
    rst.push_back(H ,"tmp.list");
    rst.push_back(G ,"tmpp.list");
#endif

#ifdef DEBUG
    Rcpp::Rcout << "phase 2" << std::endl;
#endif
    std::vector<arma::mat> gamma_new(gamma.size());
    for (int i = 0; i < gamma.size(); i++) gamma_new[i] = arma::mat(gamma[i]); // has the same structure as ga
    
    std::vector<arma::cube> V_tid(R);
    std::vector<arma::mat>  v_tid(R);
    std::vector<arma::vec>  gamma_accept(R);
#ifdef DEBUG
    std::vector<std::vector<arma::rowvec>>  gamma_props(R, std::vector<arma::rowvec>(num_cluster)); // for comparison
    std::vector<std::vector<std::vector<arma::rowvec>>>  g_props  (R, std::vector<std::vector<arma::rowvec>>(num_cluster,std::vector<arma::rowvec>(N))); // for comparison
    std::vector<std::vector<std::vector<arma::rowvec>>>  gs       (R, std::vector<std::vector<arma::rowvec>>(num_cluster,std::vector<arma::rowvec>(N))); // for comparison
    std::vector<std::vector<std::vector<arma::rowvec>>>  gt_props (R, std::vector<std::vector<arma::rowvec>>(num_cluster,std::vector<arma::rowvec>(N))); // for comparison
    std::vector<std::vector<std::vector<arma::rowvec>>>  gts      (R, std::vector<std::vector<arma::rowvec>>(num_cluster,std::vector<arma::rowvec>(N))); // for comparison
    std::vector<std::vector<std::vector<double>>>  phi_invs      (R, std::vector<std::vector<double>>(num_cluster,std::vector<double>(N))); // for comparison
    std::vector<std::vector<std::vector<double>>>  cmps          (R, std::vector<std::vector<double>>(num_cluster,std::vector<double>(N))); // for comparison
    std::vector<std::vector<std::vector<double>>>  dot1          (R, std::vector<std::vector<double>>(num_cluster,std::vector<double>(N))); // for comparison
    std::vector<std::vector<std::vector<double>>>  dot2          (R, std::vector<std::vector<double>>(num_cluster,std::vector<double>(N))); // for comparison
    std::vector<std::vector<std::vector<double>>>  sum_tmps      (R, std::vector<std::vector<double>>(num_cluster,std::vector<double>(N))); // for comparison
    std::vector<std::vector<std::vector<double>>>  sum_tmp_props (R, std::vector<std::vector<double>>(num_cluster,std::vector<double>(N))); // for comparison
    std::vector<std::vector<double>>  rrs      (R, std::vector<double>(num_cluster)); // for comparison
    std::vector<std::vector<double>>  rr_props (R, std::vector<double>(num_cluster)); // for comparison
#endif
    for (int r = 0; r < R; r++) {
      gamma_accept[r] = arma::vec(num_cluster);
             V_tid[r] = arma::cube(p[r],p[r], num_cluster);
             v_tid[r] = arma::mat(num_cluster, p[r]);
      arma::vec       y = ((Rcpp::DataFrame)dat(r))["y"];
      arma::vec       t = ((Rcpp::DataFrame)dat(r))["time"];
      arma::vec      id = ((Rcpp::DataFrame)dat(r))["id"];
      arma::uvec n_obss =             n_obs.col(r);

      for(int k = 0; k < num_cluster; k++) {
        V_tid[r].slice(k) = arma::inv(arma::inv(V0[r].slice(k)) + H[r].slice(k));
        v_tid[r].row(k)   = gamma[r].row(k) + (V_tid[r].slice(k) * (G[r].row(k).t() - arma::inv(V0[r].slice(k)) * (gamma[r].row(k) - v0[r].row(k)).t())).t();

        Rcpp::Environment MASS = Rcpp::Environment::namespace_env("MASS");
        Rcpp::Function mvrnorm = MASS["mvrnorm"];
#ifdef NORAND
        set_seed(seed_initial + iter + r+1 + k+1);
#endif
        arma::rowvec gamma_prop = Rcpp::as<arma::rowvec>(mvrnorm(
          Rcpp::_["n"]     = 1,
          Rcpp::_["mu"]    = v_tid[r].row(k),
          Rcpp::_["Sigma"] = c_gamma_tuning[r,k]*V_tid[r].slice(k)));  //  proposed new values 
#ifdef DEBUG
        gamma_props[r][k] = gamma_prop;
#endif
        double sum_tmp = 0;
        double sum_tmp_prop = 0;

        for (int i = 0; i < N; i++) {
          // TODO break unless zz_local(s,i)==j+1
          
          arma::vec  in_i = arma::conv_to<arma::vec>::from(id==i+1); // c++ index begins from 0
          arma::uvec in_i_idx = find(in_i);
          arma::vec  t_in_i = t.elem(in_i_idx);
          arma::vec one(t_in_i.size(), arma::fill::ones);
          arma::mat x_in_i = arma::join_rows(one, t_in_i, pow(t_in_i, 2), pow(t_in_i, 3));
          arma::mat y_in_i = y.elem(in_i_idx);
          arma::mat x  = x_in_i.cols(0, p[r]-1);
          arma::mat Z  = x_in_i.cols(0, q[r]-1);

          arma::rowvec gamma_prop_(gamma_prop);
          arma::rowvec gamma_(gamma[r].row(k));
          gamma_prop_.resize(p[r]);
          gamma_     .resize(p[r]);
          arma::rowvec eta_prop = gamma_prop_ * x.t() + beta[r].row(i) * Z.t();  
          arma::rowvec eta      = gamma_      * x.t() + beta[r].row(i) * Z.t();
#ifdef DEBUG
          g_props[r][k][i] = eta_prop;
          gs     [r][k][i] = eta;
#endif
          arma::rowvec q;
          arma::rowvec q_prop;
               if (dist[r].compare("gaussian") == 0)   {q = 0.5*eta%eta;      q_prop = 0.5*eta_prop%eta_prop;}
          else if (dist[r].compare("poisson" ) == 0)   {q = exp(eta);         q_prop = exp(eta_prop);}
          else if (dist[r].compare("binomial") == 0)   {q = log(1+exp(eta));  q_prop = log(1+exp(eta_prop));}
#ifdef DEBUG 
          gts     [r][k][i] = q;
          gt_props[r][k][i] = q_prop;
#endif    
          // to calculate the accept-reject function [targeted distribution]
          sum_tmp_prop += (zz_local(r,i)==k+1) * (dot(y_in_i, eta_prop) - arma::sum(q_prop));
          sum_tmp      += (zz_local(r,i)==k+1) * (dot(y_in_i, eta)      - arma::sum(q));
#ifdef DEBUG
          phi_invs     [r][k][i] = pow(phi[r][k], -1);
          cmps         [r][k][i] = zz_local(r,i)==k+1;
          dot1         [r][k][i] = dot(y_in_i, eta_prop);
          dot2         [r][k][i] = dot(y_in_i, eta);
          sum_tmps     [r][k][i] = sum_tmp;
          sum_tmp_props[r][k][i] = sum_tmp_prop;
#endif
        }

        sum_tmp_prop *= pow(phi[r][k], -1);
        sum_tmp      *= pow(phi[r][k], -1);

        double rr_prop = sum_tmp_prop - 0.5 * dot(gamma_prop     , arma::inv(V0[r].slice(k)) * gamma_prop.t());
        double rr      = sum_tmp      - 0.5 * dot(gamma[r].row(k), arma::inv(V0[r].slice(k)) * gamma[r].row(k).t());
#ifdef DEBUG
        rrs     [r][k] = rr;     
        rr_props[r][k] = rr_prop;
#endif
        double aa      = std::min(1.0, exp(rr_prop-rr));
    //     if (iter %% per == 0 & print.info == "TRUE") {
    //       cat(paste(rep('-',20),sep='',collapse=''), '\n'); 
    //       cat("dist = ", dist[[s]] , "\n")  
    //       cat(paste(rep('-',20),sep='',collapse=''), '\n'); 
    //       cat("zz.local = ", table(zz.local[[s]]), "\n")  
    //       cat("c.ga = ", c.ga[[s]], "\n")  
    //       cat("omega0.tid = ", omega0.tid[[s]], "\n")  
    //       cat("w0.tid = ", w0.tid[[s]], "\n")  
    //       cat("rr.prop = ", rr.prop, "\n")  
    //       cat("rr = ", rr, "\n")  
    //       cat("myratio = ", myratio, "\n")  
    //       cat("aa = ", aa, "\n")
    //     }
        
        if (isnan(aa)) {
          gamma_new[r].row(k) = gamma[r].row(k); gamma_accept[r].row(k) = 0;
        } else if (isinf(aa)) { // not possible since aa is at most 1.0
          gamma_new[r].row(k) = gamma_prop;      gamma_accept[r].row(k) = 1;
        } else if( arma::randu(arma::distr_param(0,1)) <= aa ) {
          gamma_new[r].row(k) = gamma_prop;      gamma_accept[r].row(k) = 1;
        } else {
          gamma_new[r].row(k) = gamma[r].row(k); gamma_accept[r].row(k) = 0;
        }
      }
    }
    gamma = gamma_new;
    // if (iter %% per == 0 & print.info == "TRUE"){print(gamma)}
#ifdef DEBUG
    rst.push_back(V_tid ,"omega0.tid");
    rst.push_back(v_tid ,"w0.tid");
    rst.push_back(gamma_props ,"ga.props");
    rst.push_back(g_props ,"g.props");
    rst.push_back(gs ,"gs");
    rst.push_back(gt_props ,"gt.props");
    rst.push_back(gts ,"gts");
    rst.push_back(phi_invs ,"phi.invs");
    rst.push_back(cmps ,"cmps");
    rst.push_back(dot1 ,"dot1");
    rst.push_back(dot2 ,"dot2");
    rst.push_back(sum_tmps ,"sum.tmps");
    rst.push_back(sum_tmp_props ,"sum.tmp.props");
    rst.push_back(gamma_accept ,"ga.accept");
    rst.push_back(gamma_new ,"ga.new");
#endif

    //-------------------------------------------------------------#
    // Adding Order Contraint 
    // Cluster with the smallest intercept will be the first group
    // Cluster with the largest intercept will be the last group
    //-------------------------------------------------------------#
    //Rcpp::Rcout << "Adding Order Contraint" << std::endl;
    //for (s in 1:R){gamma[[s]] <- gamma[[s]][order(gamma[[s]][,1]),]}

    //--------------------------------------------------------------#
    // Sample Sigma's (Random effect variances)
    //--------------------------------------------------------------#
#ifdef DEBUG
    Rcpp::Rcout << "Sample Sigma's (Random effect variances)" << std::endl;
#endif
    if (sig_var) {
      std::vector<arma::cube> tmp_list   (R);
      std::vector<arma::cube> Lambda_tid(R);
      std::vector<arma::cube> sigma_sq_u_inv(R);
      for (int r = 0; r < R; r++) {
        tmp_list[r] = arma::cube(q[r],q[r],num_cluster);
        Lambda_tid[r] = arma::cube(q[r],q[r],num_cluster);
        sigma_sq_u_inv[r] = arma::cube(q[r],q[r],num_cluster);

        for (int i = 0; i < N; i++) {
          for (int k = 0; k < num_cluster; k++) {
            tmp_list[r].slice(k) += (zz_local(r,i)==k+1) * beta[r].row(i).subvec(0,q[r]) * beta[r].row(i).subvec(0,q[r]).t();
          }
        }
        Lambda_tid[r] = lambda0[r]*Lambda0[r] + tmp_list[r];
        Rcpp::Function rWishart("rWishart");
        for (int k = 0; k < num_cluster; k++) {
          sigma_sq_u_inv[r].slice(k) = Rcpp::as<arma::mat>(rWishart(1,
            arma::mat(1,1,arma::fill::value(sum(zz_local.row(r)==k+1) + lambda0(r,k))),
            arma::inv(Lambda_tid[r].slice(k))));
          sigma_sq_u[r].slice(k) = arma::inv(sigma_sq_u_inv[r].slice(k));
        }
      }
    } else {
      for (int r = 0; r < R; r++) {
        for (int k = 0; k < num_cluster; k++) {
          arma::vec zeta_sq(q[r]);
          for (int c = 0; c < q[r]; c++) {
            arma::vec in_i = arma::conv_to<arma::vec>::from(zz_local.row(r)==k+1);
            arma::uvec in_i_idx = find(in_i);
            arma::colvec beta_ = beta[r].col(c);
            beta_ = beta_(in_i_idx);

            double c0_tid = c0(r,k) + arma::sum(in_i)/2.0;
            double d0_tid = d0(r,k) + arma::dot(beta_, beta_)/2.0;
#ifdef NORAND
            set_seed(seed_initial + iter + r+1 + k+1 + c+1);
#endif
            Rcpp::Environment MCMCpack = Rcpp::Environment::namespace_env("MCMCpack");
            Rcpp::Function rinvgamma = MCMCpack["rinvgamma"];
            zeta_sq[c] = *REAL(rinvgamma(1, c0_tid, d0_tid));
          }

          sigma_sq_u[r].slice(k) *= 0;
          sigma_sq_u[r].slice(k).diag() += zeta_sq;
        }
      }
    }
#ifdef DEBUG
    rst.push_back(sigma_sq_u ,"sigma.sq.u");
#endif

    //--------------------------------------------------------------#
    // Sample random effects via  MH algorithm 
    //--------------------------------------------------------------#
#ifdef DEBUG
    Rcpp::Rcout << "Sample random effects via  MH algorithm" << std::endl;
#endif
    std::vector<arma::cube> H2(R);
    std::vector<arma::mat>  G2(R);
    for (int r = 0; r < R; r++) {
      H2[r] = arma::cube(q[r],q[r],num_cluster);
      G2[r] = arma::mat(num_cluster,q[r]);
      arma::vec       y = ((Rcpp::DataFrame)dat(r))["y"];
      arma::vec       t = ((Rcpp::DataFrame)dat(r))["time"];
      arma::vec      id = ((Rcpp::DataFrame)dat(r))["id"];
      arma::uvec n_obss =             n_obs.col(r);

      for (int i = 0; i < N; i++) {
        arma::vec  in_i = arma::conv_to<arma::vec>::from(id==i+1); // c++ index begins from 0
        arma::uvec in_i_idx = find(in_i);
        arma::vec  t_in_i = t.elem(in_i_idx);
        arma::vec one(t_in_i.size(), arma::fill::ones);
        arma::mat x_in_i = arma::join_rows(one, t_in_i, pow(t_in_i, 2), pow(t_in_i, 3));
        arma::mat y_in_i = y.elem(in_i_idx);
        arma::mat x  = x_in_i.cols(0, p[r]-1);
        arma::mat Z  = x_in_i.cols(0, q[r]-1);

        for (int k = 0; k < num_cluster; k++) {
          arma::rowvec gamma_(gamma[r].row(k).subvec(0, p[r]-1));
          arma::rowvec eta = gamma_ * x.t() + beta[r].row(i) * Z.t();
          arma::rowvec mu;
          arma::mat kappa(x.n_rows, x.n_rows);

          if (dist[r].compare("gaussian") == 0) {
            mu = eta;
            kappa.diag() += sigma_sq_e[r][k];
          } else if (dist[r].compare("poisson") == 0) {
            mu = exp(eta); // TODO is this right?
            kappa.diag() += exp(eta);
          } else if (dist[r].compare("binomial") == 0) {
            mu = exp(eta)/(1+exp(eta)); // TODO is this right?
            arma::rowvec p = exp(eta)/pow(1+exp(eta),2); // TODO why is there a square
            kappa.diag() += p % (1-p); // element-wise product
          }

          H2[r].slice(k) +=  pow(phi[r][k],-2) * (zz_local(r,i)==k+1)*(Z.t() * kappa * Z);
          G2[r].row(k)   += (pow(phi[r][k],-1) * (zz_local(r,i)==k+1)* Z.t() * (y_in_i.t() - mu).t()).t();
        }
      }
    }
#ifdef DEBUG
    rst.push_back(H2 ,"tmp.list2");
    rst.push_back(G2 ,"tmpp.list2");
#endif

#ifdef DEBUG
    Rcpp::Rcout << "phase 2" << std::endl;
#endif
    std::vector<arma::mat> beta_new(R);
    for (int r = 0; r < R; r++) beta_new[r] = arma::mat(beta[r]); // has the same structure as beta
    
    std::vector<arma::cube> Sigma_tid  (R);
    std::vector<arma::mat>  mu_tid     (R);
    std::vector<arma::vec>  beta_accept(R);
#ifdef DEBUG
    std::vector<std::vector<std::vector<arma::rowvec>>>  beta_props(R, std::vector<std::vector<arma::rowvec>>(N, std::vector<arma::rowvec>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<arma::rowvec>>>  g_props2  (R, std::vector<std::vector<arma::rowvec>>(N,std::vector<arma::rowvec>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<arma::rowvec>>>  gs2       (R, std::vector<std::vector<arma::rowvec>>(N,std::vector<arma::rowvec>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<arma::rowvec>>>  gt_props2 (R, std::vector<std::vector<arma::rowvec>>(N,std::vector<arma::rowvec>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<arma::rowvec>>>  gts2      (R, std::vector<std::vector<arma::rowvec>>(N,std::vector<arma::rowvec>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<double>>>  phi_invs2      (R, std::vector<std::vector<double>>(N,std::vector<double>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<double>>>  cmps2          (R, std::vector<std::vector<double>>(N,std::vector<double>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<double>>>  dot12          (R, std::vector<std::vector<double>>(N,std::vector<double>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<double>>>  dot22          (R, std::vector<std::vector<double>>(N,std::vector<double>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<double>>>  sum_tmps2      (R, std::vector<std::vector<double>>(N,std::vector<double>(num_cluster))); // for comparison
    std::vector<std::vector<std::vector<double>>>  sum_tmp_props2 (R, std::vector<std::vector<double>>(N,std::vector<double>(num_cluster))); // for comparison
    std::vector<std::vector<double>>  rrs2      (R, std::vector<double>(N)); // for comparison
    std::vector<std::vector<double>>  rr_props2 (R, std::vector<double>(N)); // for comparison
#endif
    for (int r = 0; r < R; r++) {
      beta_accept[r] = arma::vec(N);
        Sigma_tid[r] = arma::cube(q[r],q[r],num_cluster);
           mu_tid[r] = arma::mat(num_cluster, q[r]);
      arma::vec       y = ((Rcpp::DataFrame)dat(r))["y"];
      arma::vec       t = ((Rcpp::DataFrame)dat(r))["time"];
      arma::vec      id = ((Rcpp::DataFrame)dat(r))["id"];
      arma::uvec n_obss =             n_obs.col(r);

      for (int i = 0; i < N; i++) {
        double sum_tmp = 0;
        double sum_tmp_prop = 0;

        int k; // TODO the R code is defo wrong
        arma::rowvec beta_prop;
        for(k = 0; k < num_cluster; k++) {
          Sigma_tid[r].slice(k) = arma::inv(arma::inv(sigma_sq_u[r].slice(k)) + H2[r].slice(k));
          mu_tid[r].row(k)      = beta[r].row(i) + (Sigma_tid[r].slice(k)  * (G2[r].row(k).t() - arma::inv(sigma_sq_u[r].slice(k)) * beta[r].row(i).t())).t();

          Rcpp::Environment MASS = Rcpp::Environment::namespace_env("MASS");
          Rcpp::Function mvrnorm = MASS["mvrnorm"];
#ifdef NORAND
          set_seed(seed_initial + iter + r+1 + k+1);
#endif
          beta_prop = Rcpp::as<arma::rowvec>(mvrnorm(
            Rcpp::_["n"]     = 1,
            Rcpp::_["mu"]    = mu_tid[r].row(k),
            Rcpp::_["Sigma"] = c_beta_tuning[r] * Sigma_tid[r].slice(k)));  //  proposed new values
#ifdef DEBUG
          beta_props[r][i][k] = beta_prop;
#endif

          arma::vec  in_i = arma::conv_to<arma::vec>::from(id==i+1); // c++ index begins from 0
          arma::uvec in_i_idx = find(in_i);
          arma::vec  t_in_i = t.elem(in_i_idx);
          arma::vec one(t_in_i.size(), arma::fill::ones);
          arma::mat x_in_i = arma::join_rows(one, t_in_i, pow(t_in_i, 2), pow(t_in_i, 3));
          arma::mat y_in_i = y.elem(in_i_idx);
          arma::mat x  = x_in_i.cols(0, p[r]-1);
          arma::mat Z  = x_in_i.cols(0, q[r]-1);

          arma::rowvec gamma_(gamma[r].row(k));
          gamma_.resize(p[r]);
          arma::rowvec eta      = gamma_ * x.t() + beta[r].row(i) * Z.t();  
          arma::rowvec eta_prop = gamma_ * x.t() + beta_prop      * Z.t();  
#ifdef DEBUG
          g_props2[r][i][k] = eta_prop;
          gs2     [r][i][k] = eta;
#endif

          arma::rowvec q;
          arma::rowvec q_prop;
               if (dist[r].compare("gaussian") == 0)   {q = 0.5*eta%eta;     q_prop = 0.5*eta_prop%eta_prop;}
          else if (dist[r].compare("poisson" ) == 0)   {q = exp(eta);        q_prop = exp(eta_prop);}
          else if (dist[r].compare("binomial") == 0)   {q = log(1+exp(eta)); q_prop = log(1+exp(eta_prop));}
#ifdef DEBUG
          gts2     [r][i][k] = q;
          gt_props2[r][i][k] = q_prop;
#endif

          sum_tmp      += pow(phi[r][k], -1) * (zz_local(r,i)==k+1) * (dot(y_in_i, eta)      - sum(q));
          sum_tmp_prop += pow(phi[r][k], -1) * (zz_local(r,i)==k+1) * (dot(y_in_i, eta_prop) - sum(q_prop));
#ifdef DEBUG
          phi_invs2     [r][i][k] = pow(phi[r][k], -1);
          cmps2         [r][i][k] = zz_local(r,i)==k+1;
          dot12         [r][i][k] = dot(y_in_i, eta_prop);
          dot22         [r][i][k] = dot(y_in_i, eta);
          sum_tmps2     [r][i][k] = sum_tmp;
          sum_tmp_props2[r][i][k] = sum_tmp_prop;
#endif
        }

        double rr_prop = sum_tmp_prop - 0.5 * dot(beta_prop     , arma::inv(sigma_sq_u[r].slice(k-1)) * beta_prop.t());
        double rr      = sum_tmp      - 0.5 * dot(beta[r].row(i), arma::inv(sigma_sq_u[r].slice(k-1)) * beta[r].row(i).t());
#ifdef DEBUG
        rrs2     [r][i] = rr;     
        rr_props2[r][i] = rr_prop;
#endif

#ifdef NORAND
        set_seed(seed_initial + iter + r+1 + i+1);
#endif
        double aa      = std::min(1.0, exp(rr_prop-rr));
        //if ( arma::randu(arma::distr_param(0,1)) < aa ) {
        Rcpp::Function runif("runif");
        if ( *REAL(runif(1,0,1)) < aa ) {
          beta_new[r].row(i) = beta_prop;      beta_accept[r].row(i) = 1;
        } else {
          beta_new[r].row(i) = beta[r].row(i); beta_accept[r].row(i) = 0;
        }
      }
    }
    beta = beta_new;
#ifdef DEBUG
    rst.push_back(Sigma_tid ,"Sigma.tid");
    rst.push_back(mu_tid ,"mu.tid");
    rst.push_back(beta_props ,"theta.props");
    rst.push_back(g_props2 ,"g.props2");
    rst.push_back(gs2 ,"gs2");
    rst.push_back(gt_props2 ,"gt.props2");
    rst.push_back(gts2 ,"gts2");
    rst.push_back(phi_invs2 ,"phi.invs2");
    rst.push_back(cmps2 ,"cmps2");
    rst.push_back(dot12 ,"dot12");
    rst.push_back(dot22 ,"dot22");
    rst.push_back(sum_tmps2 ,"sum.tmps2");
    rst.push_back(sum_tmp_props2 ,"sum.tmp.props2");
    rst.push_back(rrs2 ,"rrs2");
    rst.push_back(rr_props2 ,"rr.props2");
    rst.push_back(beta_accept ,"theta.accept");
    rst.push_back(beta_new ,"theta.new");
    rst.push_back(beta, "theta");
#endif

    //--------------------------------------------------------------#
    // Sample sigmas (residual variances)
    //--------------------------------------------------------------#
#ifdef DEBUG
    Rcpp::Rcout << "Sample sigmas (residual variances)" << std::endl;
#endif
    std::vector<arma::vec>  tmpp_list(R);
    std::vector<arma::vec>   tmp_list(R);
    for (int r = 0; r < R; r++) {
      if(dist[r].compare("gaussian")==0) {
        tmpp_list[r] = arma::vec(num_cluster);
         tmp_list[r] = arma::vec(num_cluster);
        arma::vec       y = ((Rcpp::DataFrame)dat(r))["y"];
        arma::vec       t = ((Rcpp::DataFrame)dat(r))["time"];
        arma::vec      id = ((Rcpp::DataFrame)dat(r))["id"];
        arma::uvec n_obss =             n_obs.col(r);
        for (int i = 0; i < N; i++) {
          arma::vec  in_i = arma::conv_to<arma::vec>::from(id==i+1); // c++ index begins from 0
          arma::uvec in_i_idx = find(in_i);
          arma::vec  t_in_i = t.elem(in_i_idx);
          arma::vec one(t_in_i.size(), arma::fill::ones);
          arma::mat x_in_i = arma::join_rows(one, t_in_i, pow(t_in_i, 2), pow(t_in_i, 3));
          arma::mat y_in_i = y.elem(in_i_idx);
          arma::mat x  = x_in_i.cols(0, p[r]-1);
          arma::mat Z  = x_in_i.cols(0, q[r]-1);

          for (int k = 0; k < num_cluster; k++) {
            tmpp_list[r][k] = sum(n_obss.elem(find(zz_local.row(r)==k+1)));
            arma::rowvec gamma_(gamma[r].row(k));
            gamma_.resize(p[r]);
            arma::rowvec g = gamma_ * x.t() + beta[r].row(i) * Z.t(); // this is n_i times num.cluster dimention
            tmp_list[r][k] += (zz_local(r,i)==k+1) * pow(arma::norm(y_in_i - g.t(),"fro"), 2);
          }
        }
      }
    }
#ifdef DEBUG
    rst.push_back(tmp_list ,"tmp.list3");
    rst.push_back(tmpp_list ,"tmpp.list3");
#endif

#ifdef DEBUG
    Rcpp::Rcout << "phase 2" << std::endl;
#endif
    for (int r = 0; r < R; r++) {
      if(dist[r].compare("gaussian")==0) {
        arma::mat a0_tid = a0.row(r).t() + tmpp_list[r] / 2;
        arma::mat b0_tid = b0.row(r).t() + tmp_list [r] / 2;

#ifdef NORAND
        set_seed(seed_initial + iter + r+1);
#endif
        Rcpp::Environment MCMCpack = Rcpp::Environment::namespace_env("MCMCpack");
        Rcpp::Function rinvgamma = MCMCpack["rinvgamma"];
        if (sigma_sq_e_common) {
          sigma_sq_e[r] = arma::vec(num_cluster);
          sigma_sq_e[r] += *REAL(rinvgamma(1,
            Rcpp::_["shape"]=sum(a0_tid),
            Rcpp::_["scale"]=sum(b0_tid)));
        } else {
          for (int k = 0; k < num_cluster; k++) {
            sigma_sq_e[r][k] = *REAL(rinvgamma(1,
              Rcpp::_["shape"]=a0_tid.row(k),
              Rcpp::_["scale"]=b0_tid.row(k)));
          }
        }
      }
    }
#ifdef DEBUG
    rst.push_back(sigma_sq_e ,"sigma.sq.e");
#endif

    //--------------------------------------------------------------#
    // Update Phi (Dispersion Parameters)
    //--------------------------------------------------------------#
#ifdef DEBUG
    Rcpp::Rcout << "Update Phi (Dispersion Parameters)" << std::endl;
#endif
    for (int r = 0; r < R; r++) {
      phi[r] = arma::vec(num_cluster);
      for (int k = 0; k < num_cluster; k++) {
        if (dist[r].compare("gaussian") == 0)  phi[r][k] = sigma_sq_e[r][k];
        if (dist[r].compare("poisson")  == 0)  phi[r][k] = 1;
        if (dist[r].compare("binomial") == 0)  phi[r][k] = 1;
      }
    }
#ifdef DEBUG
    rst.push_back(phi ,"phi");
#endif
    
    if (iter >= burn_in && iter % thin == 0) {
      //--------------------------------------------------------------------#
      // storing the sample;
      PPI   = rbind(PPI,ppi);
      ZZ    = rbind(ZZ,zz);
      if (num_cluster > 1) {
        T   = abind(T, ppt, Rcpp::_["along"] = 3);
      }
      ALPHA = rbind(ALPHA,alpha);

      for (int r = 0; r < R; r++) {
        if (dist[r].compare("gaussian") == 0) {
          SIGMA_SQ_E[r].row(count)   =   sigma_sq_e[r].t();
        }
           GA_ACCEPT[r].row(count)   = gamma_accept[r].t();
        THETA_ACCEPT[r].row(count)   =  beta_accept[r].t();
        THETA       [r].slice(count) =  beta       [r]; 
        
        Rcpp::Function matrix("matrix");
        if(sig_var) {
          Rcpp::Environment MCMCpack = Rcpp::Environment::namespace_env("MCMCpack");
          Rcpp::Function vech = MCMCpack["vech"];
          SIGMA_SQ_U[r] = abind(SIGMA_SQ_U[r],matrix(apply(sigma_sq_u[r],3,vech),Rcpp::_["ncol"]=num_cluster),Rcpp::_["along"] = 3);
        } else { 
          Rcpp::Function diag("diag");
          SIGMA_SQ_U[r] = abind(SIGMA_SQ_U[r],matrix(apply(sigma_sq_u[r],3,diag),Rcpp::_["ncol"]=num_cluster),Rcpp::_["along"] = 3);
        }
        
        if (num_cluster > 1) {T_LOCAL[r].slice(count) = pt.slice(r); }
        ZZ_LOCAL[r].row(count) = zz_local.row(r);
              GA[r]            = abind(GA[r], gamma[r],Rcpp::_["along"] = 3);
      }

      count = count + 1; // put at the end since c++ index is 0-based
    }

    //---------------------------------------------------------------#
    // Check acceptance rate for random effect
    //  and adjusted accordingly when necessary
    //---------------------------------------------------------------#
    if (adaptive_tuning) {
      if( count > tuning_freq && num_cluster > 1) {
        std::vector<arma::vec>    ga_acc(R);
        arma::vec theta_acc(R);
        for (int r = 0; r < R; r++) {
             ga_acc[r] = Rcpp::as<arma::vec>(apply(GA_ACCEPT[r].rows(count - tuning_freq, count), 2, mean));
          theta_acc[r] = arma::mean(beta_accept[r]);
          if (theta_acc[r] > 0.5)    c_beta_tuning[r] = std::max(0.1,c_beta_tuning[r] + 0.1);
          if (theta_acc[r] < 0.1)    c_beta_tuning[r] = std::max(0.1,c_beta_tuning[r] - 0.1);
          for (int k = 0; k < num_cluster; k++) {
            if (ga_acc[r][k] > 0.5)  c_gamma_tuning(r,k) = std::max(0.1,c_gamma_tuning(r,k) + 0.1);
            if (ga_acc[r][k] < 0.1)  c_gamma_tuning(r,k) = std::max(0.1,c_gamma_tuning(r,k) - 0.1);
          }
        }
        // if (iter %% per == 0){ print(ga.acc); print(theta.acc);}
      }
    }
#ifdef DEBUG
    rst.push_back(c_beta_tuning ,"c.theta");
    rst.push_back(c_gamma_tuning ,"c.ga");
#endif

    iter++;
    if (iter % per == 0)
      Rcpp::Rcout << "iter = " << iter << std::endl;
    if (iter == max_iter) break;
  }

  Rcpp::Rcout << "Cpp bye" << std::endl ;
  rst.push_back(PPI         , "PPI");
  rst.push_back(ZZ          , "ZZ");
  rst.push_back(T           , "T");
  rst.push_back(ALPHA       , "ALPHA");
  rst.push_back(SIGMA_SQ_E  , "SIGMA.SQ.E");
  rst.push_back(GA_ACCEPT   , "GA.ACCEPT");
  rst.push_back(THETA_ACCEPT, "THETA.ACCEPT");
  rst.push_back(THETA       , "THETA");
  rst.push_back(SIGMA_SQ_U  , "SIGMA.SQ.U");
  rst.push_back(T_LOCAL     , "T.LOCAL");
  rst.push_back(ZZ_LOCAL    , "ZZ.LOCAL");
  rst.push_back(GA          , "GA");
  rst.push_back(iter        , "iter");
  return rst;
}
