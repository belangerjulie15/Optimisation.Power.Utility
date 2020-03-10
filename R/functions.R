
####### FUNCTONS: Optimisation.Power.Utility ########################################

######## FONCTIONNEMENT #############################################################
#(1) Après chaque modification/ajout d'une fonction, to re-write the .Rd files (description files):
#library(roxygen2)
#roxygen2::roxygenise()
#
#(2)  Pour faire apparaitre la fonction dans la Help Pages: cliquer sur «Build & Reload»
#                                       OU
# Si «Build & Reload» n'est pas disponible, alors on doit passer par le package «devtools»
#library(devtools)
#build()
#install()
######################################################################################


#1) Payoff/résultat d'un call
#' Payoff of a call option
#'
#'This function return the result of a call option
#' @param X_T The value of the underlying asset
#' @param b The strike price
#' @param a The multiplicator factor (1 most of the time)
#' @param K The constant that is added to the payoff (equal to b most of the time)
#'
#' @return The result of a call option
#' @export
#'
#' @examples payoff_call(10,1,1,1)
payoff_call<-function(X_T,b,a,K){
  payoff<-a*max(X_T-b,0)+K
  return(payoff)
}

#2) Power Utility function
#'Power Utility function
#'
#' @param X_Tu Value from which we want to calculate the utility
#' @param gamma Parameter of the Power utility function. Must be different than 0
#'
#' @return The (Power) utility value corresponding to inputs
#' @export
#'
#' @examples P_Utility(10,2)
P_Utility<-function(X_Tu,gamma){
  res<-(X_Tu)^(1-gamma)/(1-gamma)
  return(res)
}

#3) Trouve x concavifié
#' Concavification Technique
#'
#' Finds the value of the concavification point in Carpenter (2000); the point that makes an optimisation function concave again
#' @param b_o The strike price of the call option (variable annuity)
#' @param a_o The multiplicator factor of the call option (variable annuity)
#' @param K_o The constant that is added to the payoff (equal to b_o most of the time)
#' @param gamma_o  Paramater Gamma in the Power Utility function. Must be different than 0
#'
#' @return The concavification point
#' @export
#'
#' @examples Find_x_theta_PU(1,1,1,2)
Find_x_theta_PU<-function(b_o,a_o,K_o,gamma_o){
  f <- function(x) payoff_call(X_T=x,b=b_o,a=a_o,K=K_o)^(-gamma_o)*x*a_o-P_Utility(payoff_call(x,b=b_o,a=a_o,K=K_o),gamma_o)+P_Utility(K_o,gamma_o)        #-Finding the fee
  result<-uniroot(f, c(b_o,b_o+100))$root
  return(result)
}

#4) Trouver le langrangien optimal (pour call bien sûr)

#' Intermediate function
#'
#' This function is an intermediate function and is needed for the function Find_lambda
#' @param r_FL Risk neural rate (striclty less than alpha_FL ), yearly composed
#' @param alpha_FL Risky rate of return, yearly composed
#' @param sigma_FL Volatility of the risky asset, yearly composed
#' @param gamma_FL Paramater Gamma in the Power Utility function. Must be different than 0
#' @param K_FL The constant that is added to the payoff (equal to b_o most of the time)
#' @param b_FL The strike price of the call option (variable annuity)
#' @param a_FL a_o The multiplicator factor of the call option (variable annuity)
#' @param T_FL Maturity, yearly composed
#' @param lambda_FL Lambda (Lagrange multiplicator)
#'
#' @return The initial asset value given a specific value of Lagrange multiplicator
#' @export
#'
#' @examples Find_lambda(r_FL=0.02,alpha_FL=0.04,sigma_FL=0.3,gamma_FL=2,K_FL=1,b_FL=1,a_FL=1,T_FL=10,lambda_FL=0.18)
Find_lambda<-function(r_FL,alpha_FL,sigma_FL,gamma_FL,K_FL,b_FL,a_FL,T_FL,lambda_FL){
  theta_FL<-(alpha_FL-r_FL)/sigma_FL
  x_tilde_FL<-Find_x_theta_PU(b_FL,a_FL,K_FL,gamma_FL)
  indi<--log(a_FL/lambda_FL*(a_FL*(x_tilde_FL-b_FL)+K_FL)^(-gamma_FL))

  tempo<-(r_FL+0.5*theta_FL^2)*T_FL
  pro_1<-pnorm((tempo-T_FL*(1-1/gamma_FL)*theta_FL^2-indi)/(theta_FL*sqrt(T_FL)))
  pro_2<-pnorm((tempo-T_FL*theta_FL^2-indi)/(theta_FL*sqrt(T_FL)))

  FL_1<-lambda_FL^(-1/gamma_FL)/(a_FL^(1-1/gamma_FL))*pro_1*exp(-tempo*(1-1/gamma_FL)+0.5*T_FL*((1-1/gamma_FL)*theta_FL)^2)
  FL_2<-(b_FL-K_FL/a_FL)*pro_2*exp(-r_FL*T_FL)

  X_0<-FL_1+FL_2
  return(X_0)
}

#'Lambda Optimal
#'
#' @param r_FLo Risk neural rate (striclty less than alpha_FL ), yearly composed
#' @param alpha_FLo Risky rate of return, yearly composed
#' @param sigma_FLo  Volatility of the risky asset, yearly composed
#' @param gamma_FLo Paramater Gamma in the Power Utility function. Must be different than 0
#' @param K_FLo The constant that is added to the payoff (equal to b_o most of the time)
#' @param b_FLo The strike price of the call option (variable annuity)
#' @param a_FLo The multiplicator factor of the call option (variable annuity)
#' @param X_0_FLo Budget value, or initial value invested in the portfolio/V.A.
#' @param T_FLo  Maturity, yearly composed
#'
#' @return The optimal Lambda (Lagrange multiplicator) given a specific intial budget value (X_0_FLo)
#' @export
#'
#' @examples lambda_optimal(r_FLo=0.02,alpha_FLo=0.04,sigma_FLo=0.2,gamma_FLo=2,K_FLo=1,b_FLo=1,a_FLo=1,X_0_FLo=1.2,T_FLo=10)
lambda_optimal<-function(r_FLo,alpha_FLo,sigma_FLo,gamma_FLo,K_FLo,b_FLo,a_FLo,X_0_FLo,T_FLo){

  f <- function(x) Find_lambda(r_FL=r_FLo,alpha_FL=alpha_FLo,sigma_FL=sigma_FLo,gamma_FL=gamma_FLo,K_FL=K_FLo,b_FL=b_FLo,a_FL=a_FLo,T_FL=T_FLo,lambda_FL = x)-X_0_FLo

  lambda_result<-uniroot(f, c(0.0000000000000001,100))$root
  return(lambda_result)
}

#5) Fonction de tradding à chaque période

# Inputs:
#invest_risq_tm1-> nombres d'unités dans actif risqué avant le rebalancement au temps t
#invest_srisq_tm1-> nombre d'unités dans actif sans riqueavant le rebalancement au temps t
#Outputs:
#nombres d'unités dans actif risqué/sans risque après le rebalancement.

#' Rebalancing & Investment
#'
#' @param S_t Risky asset value at time t (time of rebalancing)
#' @param B_t Riskless asset value at time t (time of rebalancing)
#' @param invest_risq_tm1 Number of units invested in the risky asset before rebalancing (at time t-1)
#' @param invest_srisq_tm1 Number of units invested in the riskless asset before rebalancing (at time t-1)
#' @param propor Proportion of the total investment in the risky asset
#'
#' @return A vector which contains in first positionthe number of units invest in the risky asset, and in second the riskless
#' @export
#'
#' @examples Investment_Fonct(1.4,1.2,0.3,0.7,0.5)
Investment_Fonct<-function(S_t,B_t,invest_risq_tm1,invest_srisq_tm1,propor){
  Ptf_ini<-S_t*invest_risq_tm1+invest_srisq_tm1*B_t

  investment<-rbind(propor*Ptf_ini/S_t,Ptf_ini*(1-propor)/B_t)

  return(investment)
}

#6) Inverse de la fonction Power Utility
#' Inverse Function
#'
#' @param uty Utility value
#' @param gamma_p Paramater Gamma in the Power Utility function. Must be different than 0
#'
#' @return The inverse of the Power Utility function, i.e. the value of the asset that results to the corresponding utility.
#' @export
#'
#' @examples Inverse_P_Utility(-0.4,2)
Inverse_P_Utility<-function(uty,gamma_p){
  tempo<-((1-gamma_p)*uty)^(1/(1-gamma_p))
  return(tempo)
}


#7) Fonction Power Utility avec concavification
#' Concavification of the Power Utility function
#'
#' @param X_T_con Risky asset value at time t
#' @param gamma_c Paramater Gamma in the Power Utility function. Must be different than 0.
#' @param b_c The strike price of the call option (variable annuity)
#' @param a_c The multiplicator factor of the call option (variable annuity)
#' @param K_c The constant that is added to the payoff (equal to b_o most of the time)
#'
#' @return The concavified Power Utility function
#' @export
#'
#' @examples P_Utility_con(X_T_con=1.3,gamma_c=3,b_c=1,a_c=1,K_c=1)
P_Utility_con<-function(X_T_con,gamma_c,b_c,a_c,K_c){
  x_theta<-Find_x_theta_PU(b_o=b_c,a_o=a_c,K_o=K_c,gamma_o=gamma_c)
  slope<-(P_Utility(payoff_call(x_theta,b_c,a_c,K_c),gamma_c)-P_Utility(K_c,gamma_c))/(x_theta)

  if(X_T_con<=x_theta){Chord<-P_Utility(K_c,gamma_c)+X_T_con*slope}
  else{Chord<-P_Utility(payoff_call(X_T=X_T_con,b=b_c,a=a_c,K=K_c),gamma_c)}

  return(Chord)
}


#8) Dérivée de la fonction Power Utility avec concavification
#' Derivative of the Concavified Power Utility function
#'
#' @param X_T_con Risky asset value at time t
#' @param gamma_c Paramater Gamma in the Power Utility function. Must be different than 0.
#' @param b_c The strike price of the call option (variable annuity)
#' @param a_c The multiplicator factor of the call option (variable annuity)
#' @param K_c The constant that is added to the payoff (equal to b_o most of the time)
#'
#'
#' @return The derivative of the function Power Utility from which concavification has been made
#' @export
#'
#' @examples Derivee_P_Utility_con(X_T_con=1.4,gamma_c=2,b_c=1,a_c=1,K_c=1)
Derivee_P_Utility_con<-function(X_T_con,gamma_c,b_c,a_c,K_c){
  x_theta<-Find_x_theta_PU(b_o=b_c,a_o=a_c,K_o=K_c,gamma_o=gamma_c)
  slope<-(P_Utility(payoff_call(x_theta,b_c,a_c,K_c),gamma_c)-P_Utility(K_c,gamma_c))/(x_theta)

  if(X_T_con<=x_theta){Chord<-slope}
  else{Chord<-payoff_call(X_T=X_T_con,b=b_c,a=a_c,K=K_c)^(-gamma_c)*a_c}

  return(Chord)
}

#9) Somme de deux function Power Utility
#' Sum of 2 Power Utility functions (different risk aversion)
#'
#' @param gamma1 Paramater Gamma 1 in the Power Utility function. Must be different than 0.
#' @param gamma2 Paramater Gamma 2 in the Power Utility function. Must be different than 0.
#' @param b The strike price of the call option (variable annuity)
#' @param a The multiplicator factor of the call option (variable annuity)
#' @param K The constant that is added to the payoff (equal to b_o most of the time)
#' @param x Risky asset value at time t
#'
#' @return The total utility coresponding to the asset value
#' @export
#'
#' @examples PowerUtility_2(gamma1=2,gamma2=4,b=1,a=1,K=1,x=1.5)
PowerUtility_2<-function(gamma1,gamma2,b,a,K,x){
  sumUtl<-P_Utility(x,gamma1)+P_Utility(x,gamma2)
  return(sumUtl)
}

#10) Somme de trois function Power Utility
#' Sum of 3 Power Utility functions (different risk aversion)
#'
#' @param gamma1 Paramater Gamma 1 in the Power Utility function. Must be different than 0.
#' @param gamma2 Paramater Gamma 2 in the Power Utility function. Must be different than 0.
#' @param gamma3 Paramater Gamma 3 in the Power Utility function. Must be different than 0.
#' @param b The strike price of the call option (variable annuity)
#' @param a The multiplicator factor of the call option (variable annuity)
#' @param K The constant that is added to the payoff (equal to b_o most of the time)
#' @param x Risky asset value at time t
#'
#' @return The total utility coresponding to the asset value
#' @export
#'
#' @examples PowerUtility_3(gamma1=2,gamma2=3,gamma3=5,b=1,a=1,K=1,x=1.5)
PowerUtility_3<-function(gamma1,gamma2,gamma3,b,a,K,x){
  sumUtl<-P_Utility(x,gamma1)+P_Utility(x,gamma2)+P_Utility(x,gamma3)
  return(sumUtl)
}

#11) Somme de deux function Power Utility concavifiée
#' Sum of 2 concavified Power Utility functions (different risk aversion)
#'
#' @param gamma1 Paramater Gamma 1 in the Power Utility function. Must be different than 1.
#' @param gamma2 Paramater Gamma 2 in the Power Utility function. Must be different than 1.
#' @param b The strike price of the call option (variable annuity)
#' @param a The multiplicator factor of the call option (variable annuity)
#' @param K The constant that is added to the payoff (equal to b_o most of the time)
#' @param x Risky asset value at time t
#'
#' @return The total concavified utility coresponding to the asset value
#' @export
#'
#' @examples PowerUtility_con_2(gamma1=2,gamma2=4,b=1,a=1,K=1,x=1.5)
PowerUtility_con_2<-function(gamma1,gamma2,b,a,K,x){
  sumUtl<-P_Utility_con(x, gamma1, b, a, K)+P_Utility_con(x, gamma2, b, a, K)
  return(sumUtl)
}

#12) Somme de deux function Power Utility concavifiée
#' Sum of 3 concavified Power Utility functions (different risk aversion)
#'
#' @param gamma1 Paramater Gamma 1 in the Power Utility function. Must be different than 0.
#' @param gamma2 Paramater Gamma 2 in the Power Utility function. Must be different than 0.
#' @param gamma3 Paramater Gamma 3 in the Power Utility function. Must be different than 0.
#' @param b The strike price of the call option (variable annuity)
#' @param a The multiplicator factor of the call option (variable annuity)
#' @param K The constant that is added to the payoff (equal to b_o most of the time)
#' @param x Risky asset value at time t
#'
#' @return The total concavified utility coresponding to the asset value
#' @export
#'
#' @examples PowerUtility_con_3(gamma1=2,gamma2=4,gamma3=5,b=1,a=1,K=1,x=1.6)
PowerUtility_con_3<-function(gamma1,gamma2,gamma3,b,a,K,x){
  sumUtl<-P_Utility_con(x, gamma1, b, a, K)+P_Utility_con(x, gamma2, b, a, K)++P_Utility_con(x, gamma3, b, a, K)
  return(sumUtl)
}


#13) Trouve x concavifié d'une fonction d'optimisation à plusieurs assurés
#' Concavification Technique for many insured
#'
#' Finds the value of the concavification point in Carpenter (2000); the point that makes an optimisation function concave again for a insured people group optimisation
#' @param b_o The strike price of the call option (variable annuity)
#' @param a_o The multiplicator factor of the call option (variable annuity)
#' @param K_o The constant that is added to the payoff (equal to b_o most of the time)
#' @param gamma_o  . Vector of all the Gamma of each of the insured. Each gamma must be different than 0

#'
#' @return The concavification point of the optimization function
#' @export
#'
#' @examples Find_x_theta_PU_Numerous(b_o=1,a_o=1,K_o=1,vector_gamma=c(2,3,5))
Find_x_theta_PU_Numerous<-function(b_o,a_o,K_o,vector_gamma){

  f <- function(x) sum(P_Utility(payoff_call(x,b_o,a_o,K_o),vector_gamma))-sum(P_Utility(K_o,vector_gamma))-x*sum((x)^(-vector_gamma))
  sln_x<-uniroot(f, c(b_o,b_o+100))$root

  return(sln_x)
}

# verif_Find_theta<-function(x){
#   res<--x^(-1)-0.5*x^(-2)-0.25*x^(-4)-(-1-0.5-0.25)-x*(x^(-2)+x^(-3)+x^(-5))
#   return(res)
# }
# uniroot(verif_Find_theta, c(1,1+100))$root

#14) Test: fonction d'une droite
#' Droite
#'
#' @param pente Pente de la droite à tracer
#' @param origine Origine de la droite à tracer
#' @param val_x Valeur de l'axe des x
#'
#' @return Valeur sur l'axe des ordonnées correspondate
#' @export
#'
#' @examples Droite(-3, 0.2,40)
Droite<-function(pente, origine,val_x){
  return(origine+val_x*pente)
}

#15) Fonction CRRA
#' CRRA (Constant Relative Risk Aversion)
#'
#' @param x Risky asset value
#' @param gamma Paramater Gamma in the Power Utility function. Must be different than 1.
#'
#' @return Utility of a CRRA utility function
#' @export
#'
#' @examples CRRA(1.5,4)
CRRA<-function(x,gamma){
  result<-(x^(1-gamma)-1)/(1-gamma)
  return(result)
}



#16) Dérivée de la fonction Power Utility
#' Derivative of the Concavified Power Utility function
#'
#' @param x Risky asset value
#' @param gamma Paramater Gamma in the Power Utility function. Must be different than 1.
#'
#' @return Derivative of Power Utility evaluated at x.
#' @export
#'
#' @examples derivee_Power_Utility(1.5,4)
derivee_Power_Utility<-function(x,gamma){
  inv<-x^(-gamma)
  return(inv)
}
################################################################


