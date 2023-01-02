#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>


using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List sim_beta(const arma::mat obs
                    ,const int model
                    ,const int L
                    ,const int people
                    ,const int itemnum
                    ,const int testlet
                    ,const int total_items
                    ,const arma::cube alpha
                    ,const arma::cube theta
                    ,const arma::mat delta
                    ,const arma::mat tau
                    ,const arma::mat sigma_item_cor
                    ){
  // arma::cube M(L,people,total_items);
  // arma::cube N(L,people,total_items);

  arma::cube nu_obs2(L,people,total_items);
  arma::cube nu_rep2(L,people,total_items);
  double r_rep;
  double mu;
  arma::cube var(L,people,total_items);
  
  vec r_obs_minus_mu(total_items);
  vec r_rep_minus_mu(total_items);

  vec sgddm_obs = zeros<vec>(L);
  vec sgddm_rep = zeros<vec>(L);
  
  int D_1 = itemnum-1;
	
	if (model == 1){
	  
	  double M;
	  double N;
	  
	  for(int r = 0; r < L; ++r) {
	    
	    mat cov_obs = zeros<mat>(total_items,total_items);
	    mat cov_rep = zeros<mat>(total_items,total_items);
	    
	    for(int n = 0; n < people; ++n) {
	      for(int k = 0; k < total_items; ++k) {
	        M = as_scalar(exp(((sum(alpha.subcube(r, k, 0, r, k, D_1) % theta.subcube(r, n, 0, r, n, D_1), 2)-delta(r,k))+tau(r,k))/2));
	        N = as_scalar(exp((-(sum(alpha.subcube(r, k, 0, r, k, D_1) % theta.subcube(r, n, 0, r, n, D_1), 2)-delta(r,k))+tau(r,k))/2));
	        r_rep = R::rbeta(M,N);
	        
	        mu = M/(M+N);
	        var(r,n,k) = (mu*(1-mu))/(M+N+1);
	        
	        r_obs_minus_mu(k) = obs(n,k) - mu;
	        r_rep_minus_mu(k) = r_rep  - mu;
	        
	        nu_obs2(r,n,k)  = pow(r_obs_minus_mu(k), 2);
	        nu_rep2(r,n,k)  = pow(r_rep_minus_mu(k), 2);
	      }
	      
	      // sum over people
	      for(int j = 0; j < total_items; ++j) {
	        for(int i = j; i < total_items; ++i) {
	          cov_obs(i,j) = cov_obs(i,j) + r_obs_minus_mu(i) * r_obs_minus_mu(j);
	          cov_rep(i,j) = cov_rep(i,j) + r_rep_minus_mu(i) * r_rep_minus_mu(j);
	          
	          cov_obs(j,i) = cov_obs(i,j);
	          cov_rep(j,i) = cov_rep(i,j);
	        }
	      }
	      
	    } // end of people
	    
	    // diagonals are ones; off-diagonals are not correlations because mu is not same over people.
	    for(int j = 0; j < total_items; ++j) {
	      for(int i = (1+j); i < total_items; ++i) {
	        cov_obs(i,j) = cov_obs(i,j)/sqrt(cov_obs(i,i) * cov_obs(j,j));
	        cov_rep(i,j) = cov_rep(i,j)/sqrt(cov_rep(i,i) * cov_rep(j,j));
	        
	        cov_obs(j,i) = cov_obs(i,j);
	        cov_rep(j,i) = cov_rep(i,j);
	      }
	    }
	    // diagonals
	    for(int j = 0; j < total_items; ++j) {
	      cov_obs(j,j) = 1;
	      cov_rep(j,j) = 1;
	    }
	    
	    for(int j = 0; j < total_items; ++j) {
	      for(int i = (1+j); i < total_items; ++i) {
	        sgddm_obs(r) = sgddm_obs(r) + abs(cov_obs(i,j));
	        sgddm_rep(r) = sgddm_rep(r) + abs(cov_rep(i,j));
	      }
	    }
	    sgddm_obs(r) = sgddm_obs(r)*2/(total_items*(total_items-1));
	    sgddm_rep(r) = sgddm_rep(r)*2/(total_items*(total_items-1));
	    
	  } // end of replications
	  
	  return Rcpp::List::create(Rcpp::Named("nu_obs2",nu_obs2)
                               ,Rcpp::Named("nu_rep2",nu_rep2)
                               ,Rcpp::Named("var",var)
                               ,Rcpp::Named("sgddm_obs",sgddm_obs)
                               ,Rcpp::Named("sgddm_rep",sgddm_rep)
	  );
	  
	}else if (model == 2){
	  
	  vec z(itemnum);
	  mat covar(itemnum,itemnum);
	  vec SD(itemnum);
	  vec value(itemnum);
	  double itemnum2 = itemnum*itemnum;
	  
	  
	  double index;
	  
	  vec M(total_items);
	  vec N(total_items);
	  
	  vec MU = zeros<vec>(itemnum);
	  
	  for(int r = 0; r < L; ++r) {
	    
	    mat cov_obs = zeros<mat>(total_items,total_items);
	    mat cov_rep = zeros<mat>(total_items,total_items);
	    
	    for(int n = 0; n < people; ++n) {
	      for(int k = 0; k < total_items; ++k) {
	        M(k) = as_scalar(exp(((sum(alpha.subcube(r, k, 0, r, k, D_1) % theta.subcube(r, n, 0, r, n, D_1), 2)-delta(r,k))+tau(r,k))/2));
	        N(k) = as_scalar(exp((-(sum(alpha.subcube(r, k, 0, r, k, D_1) % theta.subcube(r, n, 0, r, n, D_1), 2)-delta(r,k))+tau(r,k))/2));
	      }
	      

	     

	      for(int k = 0; k < testlet; ++k) {
	        for(int j = 0; j < itemnum; ++j) {
	          for(int i = 0; i < itemnum; ++i) {
	            covar(i,j) = sigma_item_cor(r, k*itemnum2 + j*itemnum + i );
	          }
	        }
	        z = mvnrnd(MU, covar, 1);
	        SD = sqrt(covar.diag());
	        for(int i = 0; i < itemnum; ++i) {
	          value(i) = R::pnorm( z(i) , 0, SD(i), true, false);
	          index = k*itemnum + i;
	          r_rep = R::qbeta( value(i), M(index), N(index), true, false);
	          
	          mu = M(index)/(M(index) + N(index));
	          var(r,n,index) = (mu*(1-mu))/(M(index)+N(index)+1);
	          
	          r_obs_minus_mu(index) = obs(n,index) - mu;
	          r_rep_minus_mu(index) = r_rep - mu;
	          
	          nu_obs2(r,n,index) = pow(r_obs_minus_mu(index), 2);
	          nu_rep2(r,n,index) = pow(r_rep_minus_mu(index), 2);
	        }
	      } // end of testlet

        
        // sum over people
        for(int j = 0; j < total_items; ++j) {
          for(int i = j; i < total_items; ++i) {
            cov_obs(i,j) = cov_obs(i,j) + r_obs_minus_mu(i) * r_obs_minus_mu(j);
            cov_rep(i,j) = cov_rep(i,j) + r_rep_minus_mu(i) * r_rep_minus_mu(j);
            
            cov_obs(j,i) = cov_obs(i,j);
            cov_rep(j,i) = cov_rep(i,j);
          }
        }
        


	    } // end of people
	    
	    
      // diagonals are ones; off-diagonals are not correlations because mu is not same over people.
      for(int j = 0; j < total_items; ++j) {
        for(int i = (1+j); i < total_items; ++i) {
          cov_obs(i,j) = cov_obs(i,j)/sqrt(cov_obs(i,i) * cov_obs(j,j));
          cov_rep(i,j) = cov_rep(i,j)/sqrt(cov_rep(i,i) * cov_rep(j,j));
          
          cov_obs(j,i) = cov_obs(i,j);
          cov_rep(j,i) = cov_rep(i,j);
        }
      }
      // diagonals
      for(int j = 0; j < total_items; ++j) {
        cov_obs(j,j) = 1;
        cov_rep(j,j) = 1;
      }
      
      for(int j = 0; j < total_items; ++j) {
        for(int i = (1+j); i < total_items; ++i) {
          sgddm_obs(r) = sgddm_obs(r) + abs(cov_obs(i,j));
          sgddm_rep(r) = sgddm_rep(r) + abs(cov_rep(i,j));
        }
      }
      sgddm_obs(r) = sgddm_obs(r)*2/(total_items*(total_items-1));
	    sgddm_rep(r) = sgddm_rep(r)*2/(total_items*(total_items-1));
      
	    
	  } // end of replications
	  
	  return Rcpp::List::create(Rcpp::Named("nu_obs2",nu_obs2)
                               ,Rcpp::Named("nu_rep2",nu_rep2)
                               ,Rcpp::Named("var",var)
                               ,Rcpp::Named("sgddm_obs",sgddm_obs)
                               ,Rcpp::Named("sgddm_rep",sgddm_rep)
	  );
	}
	
	


	
}

