# include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
#include <math.h>
#include <numeric>
#include "c_which.h"

using namespace Rcpp;
using namespace sugar;
using namespace arma;

//arma::vec  c_which(int input_id, arma::vec id,  arma::vec n_obs, NumericVector indata);

// [[Rcpp::export()]]
arma:: mat  LL(const Rcpp::List &fit, int fast_version) {   // calculate likelihood based on the BCC object

  //Rcpp::Environment myEnv = Rcpp::Environment::global_env();
  //Function c_which = myEnv["c_which"];
  Function table("table");
  Function prod("prod");
  //Rcpp::Environment base("package:base");
    //Function table = base["table"];
    //Function prod = base["prod"];
    //Function is_finite = base["is.finite"];

  Function dpois("dpois");
  Function dbinom("dbinom");
  //Rcpp::Environment stats("package:stats");
   //Function dpois = stats["dpois"];
   //Function dbinom = stats["dbinom"];

   //Function dmvnorm("dmvnorm");
   //Function mvrnorm("mvrnorm");

Rcpp::Environment mclust = Rcpp::Environment::namespace_env("mclust");
Rcpp::Function dmvnorm = mclust["dmvnorm"];

    //Function dmvnorm = mclust["dmvnorm"];
Rcpp::Environment MASS = Rcpp::Environment::namespace_env("MASS");
Rcpp::Function mvrnorm = MASS["mvrnorm"];
    //Function mvrnorm = MASS["mvrnorm"];

  int max_iter = fit["max.iter"];
  int burn_in = fit["burn.in"];
  int thin = fit["thin"];

  int N = fit["N"];
  int R = fit["R"];
  Rcpp:: List K = fit["K"];
  Rcpp:: List k = fit["k"];

  CharacterVector  dist = fit["dist"];

  int num_cluster = fit["num.cluster"];

  Rcpp:: List  dat = fit["dat"];
  Rcpp:: List  summary_stat = fit["summary.stat"];
  Rcpp:: List  ZZ_LOCAL = fit["ZZ.LOCAL"];
  Rcpp:: List  GA = fit["GA"];
  Rcpp:: List  SIGMA_SQ_E = fit["SIGMA.SQ.E"];
  Rcpp:: List  SIGMA_SQ_U = fit["SIGMA.SQ.U"];


  Rcpp:: List  THETA = fit["THETA"];
  arma:: mat ppi = zeros(1,num_cluster);

  int num_sample  =  (max_iter - burn_in)/thin;
  if (fast_version==1) {num_sample = std::min(num_sample,100); }

  arma:: mat  LOG_LIK_ITER = zeros(N,num_sample);

  for (int count = 0; count < num_sample; ++count){  // Loop MCMC

    if (num_cluster==1)  {
      ppi(0,0) = 1;
      }
    else{
      arma:: mat  PPI = fit["PPI"]; ppi.row(0) = PPI.row(count);
      }

    arma:: vec log_lik_iter(N);log_lik_iter.fill(0);


   for (int i = 0; i < N; ++i) { // Loop individual

     double tmp = 0;

     for (int j = 0; j < num_cluster; ++j) { // Loop cluster

       arma:: vec ffs(R);ffs.fill(1);
       arma:: vec vv(R); vv.fill(1);

       double ff = 1;

       for (int s = 0; s < R; ++s) {  // Loop marker
         arma::mat  ZZ_LOCAL_s = ZZ_LOCAL(s);
         arma::mat  ZZ_LOCAL_s_count = ZZ_LOCAL_s.row(count);
         arma::cube  GA_s = GA(s);
         arma::mat GA_s_count = GA_s.slice(count);
         arma::mat  SIGMA_SQ_E_s = SIGMA_SQ_E(s);
         arma::mat SIGMA_SQ_E_s_count = SIGMA_SQ_E_s.row(count);
         arma::cube  SIGMA_SQ_U_s = SIGMA_SQ_U(s);
         arma::mat SIGMA_SQ_U_s_count = SIGMA_SQ_U_s.slice(count);
         arma::cube  THETA_s = THETA(s);
         arma::mat THETA_s_count = THETA_s.slice(count);
         Rcpp::List mydat = dat(s);
         arma:: vec y_s = mydat["y"];
         arma:: vec xt = mydat["time"];
         arma:: vec ids = mydat["id"];
         arma:: vec n_obss = as<arma::vec>(table(ids)); // obtain number of observations

         arma::vec  xt_i = c_which(i+1,ids,n_obss,Rcpp::NumericVector(Rcpp::wrap(xt)));

         arma::vec  y_si = c_which(i+1,ids,n_obss,Rcpp::NumericVector(Rcpp::wrap(y_s)));


         arma::vec v_one_i (n_obss(i));  v_one_i.fill(1);

         arma::mat m = join_horiz(join_horiz(v_one_i,
                                             join_horiz(xt_i,pow(xt_i,2))),
                                             pow(xt_i,3)) ;
         arma::mat mz = join_horiz(join_horiz(v_one_i,
                                              join_horiz(xt_i,pow(xt_i,2))),
                                              pow(xt_i,3)) ;

         arma::mat  g =  m.submat(0,0,n_obss(i)-1,int(k(s))-1)*
           trans(GA_s_count.submat(j,0,j,int(k(s))-1)) +
           m.submat(0,0,n_obss(i)-1,int(K(s))-1)*
           trans(THETA_s_count.submat(i,0,i,int(K(s))-1));

         if (dist(s)=="gaussian"){
           arma::vec  sigma_sq_e_s_l(n_obss(i));
           sigma_sq_e_s_l.fill(SIGMA_SQ_E_s_count(j));
           ff = as<double>(dmvnorm(_("data")= trans(y_si),
                                   _("mean") = trans(g),
                                   _("sigma")= diagmat(sigma_sq_e_s_l)));
           }
         if (dist(s)=="poisson"){
           ff = as<double>(prod(as<arma::vec>(dpois(y_si,exp(g)))));
           }
         if (dist(s)=="binomial"){
           arma::mat mu = exp(g)/(1+exp(g));
           ff = as<double>(prod(as<arma::vec>(dbinom(y_si,_("size")=1,
                                                     _("prob")=mu))));
           }
         ffs(s) = ff;

         if (num_cluster==1) {
           vv(s) = 1;
           }
         else{
           arma:: vec alpha = fit["alpha"];
           vv(s) = (ZZ_LOCAL_s_count(0,i)==(j+1))*alpha(s) +
             (ZZ_LOCAL_s_count(0,i)!=(j+1))*(1-alpha(s))/(num_cluster-1);
           }
       }  // Closed Loop marker
     double tmp_prod1 = as<double>(prod(vv));
     double tmp_prod2 = as<double>(prod(ffs));
     tmp  =  tmp + (tmp_prod1*ppi(0,j))*tmp_prod2;

   }// Closed Loop cluster

  // CharacterVector ev = as<CharacterVector>(is_finite(log(tmp)));
  // if (ev(0)=="TRUE") {log_lik_iter(i) = log(tmp);} else{log_lik_iter(i) = 0;}
  log_lik_iter(i) = log(tmp);

  } // Closed Loop Individual
 LOG_LIK_ITER.col(count) = log_lik_iter;

 }// Closed Loop MCMC

 return LOG_LIK_ITER;
}












