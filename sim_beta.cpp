#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>


using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List Cox_Snell_fit_rcpp(
                     const arma::vec t
                    ,const arma::cube r_Cox_Snell
                    ,const arma::mat eap_r_Cox_Snell
                    ){
	int L = r_Cox_Snell.n_rows;
	int people = r_Cox_Snell.n_cols;
	int total_items = r_Cox_Snell.n_slices;
	
	int t_L = t.n_elem;
	arma::cube lambda_item(t_L,L,total_items);
	arma::cube lambda_person(t_L,L,people);
	
	arma::mat lambda_item_eap(t_L,total_items);
	arma::mat lambda_person_eap(t_L,people);
	
	
	for(int k = 0; k < t_L; ++k) {	
		for(int r = 0; r < L; ++r) {
			//
			for(int n = 0; n < people; ++n) {
				for(int i = 0; i < total_items; ++i) {	                 
                    lambda_person(k,r,n) += r_Cox_Snell(r,n,i) <= t(k);			
				}
				lambda_person(k,r,n) /= total_items;
				if (lambda_person(k,r,n)==1.0){
					lambda_person(k,r,n) = NA_REAL;
				}else{
					lambda_person(k,r,n) = -log(1 - lambda_person(k,r,n));
				}				
			}
			//
			for(int i = 0; i < total_items; ++i) {
				for(int n = 0; n < people; ++n) {
					lambda_item(k,r,i) += r_Cox_Snell(r,n,i) <= t(k);	
				}
				lambda_item(k,r,i) /= people;				
				if (lambda_item(k,r,i)==1.0){
					lambda_item(k,r,i) = NA_REAL;
				}else{
					lambda_item(k,r,i) = -log(1 - lambda_item(k,r,i));
				}
			}
			
			
		}
		//
		for(int n = 0; n < people; ++n) {
			for(int i = 0; i < total_items; ++i) {	                 
				lambda_person_eap(k,n) += eap_r_Cox_Snell(n,i) <= t(k);			
			}
			lambda_person_eap(k,n) /= total_items;			
			if (lambda_person_eap(k,n)==1.0){
				lambda_person_eap(k,n) = NA_REAL;
			}else{
				lambda_person_eap(k,n) = -log(1 - lambda_person_eap(k,n));
			}
		}		
		//
		for(int i = 0; i < total_items; ++i) {
			for(int n = 0; n < people; ++n) {
				lambda_item_eap(k,i) += eap_r_Cox_Snell(n,i) <= t(k);	
			}
			lambda_item_eap(k,i) /= people;			
			if (lambda_item_eap(k,i)==1.0){
				lambda_item_eap(k,i) = NA_REAL;
			}else{
				lambda_item_eap(k,i) = -log(1 - lambda_item_eap(k,i));
			}
		}		
	}

    return Rcpp::List::create(  Rcpp::Named("t",t)
	                           ,Rcpp::Named("lambda_item",lambda_item)
                               ,Rcpp::Named("lambda_person",lambda_person)
                               ,Rcpp::Named("lambda_item_eap",lambda_item_eap)
                               ,Rcpp::Named("lambda_person_eap",lambda_person_eap)	
	);
	

}		

// [[Rcpp::export]]
double plogitnormal_rcpp(const double x
                        ,const double mu
						,const double var
						 ){
	double logit_x = log(x/(1-x));
	double v = (logit_x-mu)/sqrt(2*var);
	double v2 = pow(v, 2);	
	double erf_v = R::pchisq( 2 * v2, 1, true, false ) * sign(v); // same with pracma::erf(v)
	double u = 0.5*(1 + erf_v);
	return(u);
}

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
  arma::cube prob(L,people,total_items);
  arma::cube r_rep(L,people,total_items);
  double mu;
  arma::cube var(L,people,total_items);
  
  vec r_obs_minus_mu(total_items);
  vec r_rep_minus_mu(total_items);

  vec sgddm_obs = zeros<vec>(L);
  vec sgddm_rep = zeros<vec>(L);
  
  arma::cube r_Cox_Snell(L,people,total_items);
  
  int loc = itemnum-1;

  cube pred_lower_yi = zeros<cube>(L,people,total_items);
  
	
	if (model == 0){
		
		int cnt = 0;
		int cnt2 = 0;
        vec VAR(itemnum);
		mat covar(itemnum,itemnum);
		double itemnum2 = itemnum*itemnum;

		double u;
		
		double x;
		double logit_x;
		
		vec z(itemnum);
		vec mu = zeros<vec>(itemnum);
		
		for(int r = 0; r < L; ++r) {
			for(int n = 0; n < people; ++n) {
				cnt = 0;
				cnt2 = 0;
				for(int k = 0; k < testlet; ++k) {
					for(int j = 0; j < itemnum; ++j) {
					  for(int i = 0; i < itemnum; ++i) {
						covar(i,j) = sigma_item_cor(r, k*itemnum2 + j*itemnum + i );
					  }
					}					
					VAR = covar.diag(0);				
					for(int i = 0; i < itemnum; ++i) {					
						mu(i) = delta(r,cnt) + as_scalar(sum(alpha.subcube(r, cnt, 0, r, cnt, loc) % theta.subcube(r, n, 0, r, n, loc), 2));
						u = plogitnormal_rcpp(obs(n,cnt),mu(i),VAR(i));
						
						//for Cox-Snell plots
						r_Cox_Snell(r,n,cnt) = -log(1 - u);
						
						// prob
						x = obs(n,cnt);
						logit_x = log(x/(1-x));
						prob(r,n,cnt) = exp( - log(x*(1-x))  + R::dnorm( logit_x, mu(i), sqrt(VAR(i)), true ) );

						cnt += 1;
					}
					
					z = mvnrnd(mu, covar, 1);
					for(int i = 0; i < itemnum; ++i) {	
						r_rep(r,n,cnt2) = 1/(1 + exp(-z(i)));
						cnt2 += 1;
					}
					
				}
			}
		}
		
	  return Rcpp::List::create(Rcpp::Named("prob",prob)
	                           ,Rcpp::Named("r_rep",r_rep)
	                           ,Rcpp::Named("nu_obs2",nu_obs2)
                               ,Rcpp::Named("nu_rep2",nu_rep2)
                               ,Rcpp::Named("var",var)
                               ,Rcpp::Named("sgddm_obs",sgddm_obs)
                               ,Rcpp::Named("sgddm_rep",sgddm_rep)
							   ,Rcpp::Named("r_Cox_Snell",r_Cox_Snell)
	  );
	  
	  
	}else if (model == 1){
	  
	  double M;
	  double N;
	  
	  cube cov_obs = zeros<cube>(L,total_items,total_items);
	  cube cov_rep = zeros<cube>(L,total_items,total_items);
	  
	  for(int r = 0; r < L; ++r) {
	    
	    for(int n = 0; n < people; ++n) {
	      for(int k = 0; k < total_items; ++k) {
	        M = as_scalar(exp(((sum(alpha.subcube(r, k, 0, r, k, loc) % theta.subcube(r, n, 0, r, n, loc), 2)-delta(r,k))+tau(r,k))/2));
	        N = as_scalar(exp((-(sum(alpha.subcube(r, k, 0, r, k, loc) % theta.subcube(r, n, 0, r, n, loc), 2)-delta(r,k))+tau(r,k))/2));
	        r_rep(r,n,k) = R::rbeta(M,N);
			prob(r,n,k) = R::dbeta(obs(n,k),M,N,false);
	        
	        mu = M/(M+N);
	        var(r,n,k) = (mu*(1-mu))/(M+N+1);
	        
	        r_obs_minus_mu(k) = obs(n,k) - mu;
	        r_rep_minus_mu(k) = r_rep(r,n,k)  - mu;
	        
	        nu_obs2(r,n,k)  = pow(r_obs_minus_mu(k), 2);
	        nu_rep2(r,n,k)  = pow(r_rep_minus_mu(k), 2);
			
			r_Cox_Snell(r,n,k) = -log(1-R::pbeta( obs(n,k), M, N, true, false));
	      }
	      
	      // sum over people
	      for(int j = 0; j < total_items; ++j) {
	        for(int i = j; i < total_items; ++i) {
	          cov_obs(r,i,j) = cov_obs(r,i,j) + r_obs_minus_mu(i) * r_obs_minus_mu(j);
	          cov_rep(r,i,j) = cov_rep(r,i,j) + r_rep_minus_mu(i) * r_rep_minus_mu(j);
	          
	          cov_obs(r,j,i) = cov_obs(r,i,j);
	          cov_rep(r,j,i) = cov_rep(r,i,j);
	        }
	      }
	      
	    } // end of people
	    
	    // diagonals are ones; off-diagonals are not correlations because mu is not same over people.
	    for(int j = 0; j < total_items; ++j) {
	      for(int i = (1+j); i < total_items; ++i) {
	        cov_obs(r,i,j) = cov_obs(r,i,j)/sqrt(cov_obs(r,i,i) * cov_obs(r,j,j));
	        cov_rep(r,i,j) = cov_rep(r,i,j)/sqrt(cov_rep(r,i,i) * cov_rep(r,j,j));
	        
	        cov_obs(r,j,i) = cov_obs(r,i,j);
	        cov_rep(r,j,i) = cov_rep(r,i,j);
	      }
	    }
	    // diagonals
	    for(int j = 0; j < total_items; ++j) {
	      cov_obs(r,j,j) = 1;
	      cov_rep(r,j,j) = 1;
	    }
	    
	    for(int j = 0; j < total_items; ++j) {
	      for(int i = (1+j); i < total_items; ++i) {
	        sgddm_obs(r) = sgddm_obs(r) + abs(cov_obs(r,i,j));
	        sgddm_rep(r) = sgddm_rep(r) + abs(cov_rep(r,i,j));
	      }
	    }
	    sgddm_obs(r) = sgddm_obs(r)*2/(total_items*(total_items-1));
	    sgddm_rep(r) = sgddm_rep(r)*2/(total_items*(total_items-1));
	    
	  } // end of replications
	  
	  return Rcpp::List::create(Rcpp::Named("prob",prob)
	                           ,Rcpp::Named("r_rep",r_rep)
	                           ,Rcpp::Named("nu_obs2",nu_obs2)
                               ,Rcpp::Named("nu_rep2",nu_rep2)
                               ,Rcpp::Named("var",var)
							   ,Rcpp::Named("cov_obs",cov_obs)
							   ,Rcpp::Named("sgddm_obs",sgddm_obs)
                               ,Rcpp::Named("cov_rep",cov_rep)
                               ,Rcpp::Named("sgddm_rep",sgddm_rep)
							   ,Rcpp::Named("r_Cox_Snell",r_Cox_Snell)
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
	  
	  vec M_eap(total_items);
	  vec N_eap(total_items);
	  
	  vec MU = zeros<vec>(itemnum);
	  
	  cube cov_obs = zeros<cube>(L,total_items,total_items);
	  cube cov_rep = zeros<cube>(L,total_items,total_items); 
	  
	  
	  for(int r = 0; r < L; ++r) {
	    

	    
	    for(int n = 0; n < people; ++n) {
	      for(int k = 0; k < total_items; ++k) {
	        M(k) = as_scalar(exp(((sum(alpha.subcube(r, k, 0, r, k, loc) % theta.subcube(r, n, 0, r, n, loc), 2)-delta(r,k))+tau(r,k))/2));
	        N(k) = as_scalar(exp((-(sum(alpha.subcube(r, k, 0, r, k, loc) % theta.subcube(r, n, 0, r, n, loc), 2)-delta(r,k))+tau(r,k))/2));
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
	          r_rep(r,n,index) = R::qbeta( value(i), M(index), N(index), true, false);
	          prob(r,n,index) = R::dbeta(obs(n,index),M(index), N(index), false );
			  
	          mu = M(index)/(M(index) + N(index));
	          var(r,n,index) = (mu*(1-mu))/(M(index)+N(index)+1);
	          
	          r_obs_minus_mu(index) = obs(n,index) - mu;
	          r_rep_minus_mu(index) = r_rep(r,n,index) - mu;
	          
	          nu_obs2(r,n,index) = pow(r_obs_minus_mu(index), 2);
	          nu_rep2(r,n,index) = pow(r_rep_minus_mu(index), 2);
			  
			  r_Cox_Snell(r,n,index) = -log(1-R::pbeta( obs(n,index), M(index), N(index), true, false));
			  
			  pred_lower_yi(r,n,index) = obs(n,index) > r_rep(r,n,index);
	        }
	      } // end of testlet

        
        // sum over people
        for(int j = 0; j < total_items; ++j) {
          for(int i = j; i < total_items; ++i) {
            cov_obs(r,i,j) = cov_obs(r,i,j) + r_obs_minus_mu(i) * r_obs_minus_mu(j);
            cov_rep(r,i,j) = cov_rep(r,i,j) + r_rep_minus_mu(i) * r_rep_minus_mu(j);
            
            cov_obs(r,j,i) = cov_obs(r,i,j);
            cov_rep(r,j,i) = cov_rep(r,i,j);
          }
        }
        


	    } // end of people
	    
	    
      // diagonals are ones; off-diagonals are not correlations because mu is not same over people.
      for(int j = 0; j < total_items; ++j) {
        for(int i = (1+j); i < total_items; ++i) {
          cov_obs(r,i,j) = cov_obs(r,i,j)/sqrt(cov_obs(r,i,i) * cov_obs(r,j,j));
          cov_rep(r,i,j) = cov_rep(r,i,j)/sqrt(cov_rep(r,i,i) * cov_rep(r,j,j));
          
          cov_obs(r,j,i) = cov_obs(r,i,j);
          cov_rep(r,j,i) = cov_rep(r,i,j);
        }
      }
      // diagonals
      for(int j = 0; j < total_items; ++j) {
        cov_obs(r,j,j) = 1;
        cov_rep(r,j,j) = 1;
      }
      
      for(int j = 0; j < total_items; ++j) {
        for(int i = (1+j); i < total_items; ++i) {
          sgddm_obs(r) = sgddm_obs(r) + abs(cov_obs(r,i,j));
          sgddm_rep(r) = sgddm_rep(r) + abs(cov_rep(r,i,j));
        }
      }
      sgddm_obs(r) = sgddm_obs(r)*2/(total_items*(total_items-1));
	    sgddm_rep(r) = sgddm_rep(r)*2/(total_items*(total_items-1));
      
	    
	  } // end of replications
	  
	  return Rcpp::List::create(Rcpp::Named("prob",prob)
	                           ,Rcpp::Named("r_rep",r_rep)
	                           ,Rcpp::Named("nu_obs2",nu_obs2)
                               ,Rcpp::Named("nu_rep2",nu_rep2)
                               ,Rcpp::Named("var",var)
							   ,Rcpp::Named("cov_obs",cov_obs)
							   ,Rcpp::Named("cov_rep",cov_rep)
                               ,Rcpp::Named("sgddm_obs",sgddm_obs)
                               ,Rcpp::Named("sgddm_rep",sgddm_rep)
							   ,Rcpp::Named("r_Cox_Snell",r_Cox_Snell)
							   ,Rcpp::Named("pred_lower_yi",pred_lower_yi)
	  );
	}
	
	


	
}

