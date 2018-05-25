//saralamba@gmail.com


#include <Rcpp.h>
#include <math.h>  
using namespace Rcpp;


// [[Rcpp::export]]
List SIRSmodel(double t, NumericVector state, NumericVector parameters) {

  //define parameters
  double R0 = parameters["R0"];
  double durinf = parameters["durinf"]; 
  double durimm = parameters["durimm"];
  double P = parameters["P"];
  double startvac = parameters["startvac"];
  double durvac = parameters["durvac"];
  double covvac = parameters["covvac"];
  double starttreat = parameters["starttreat"];
  double covtreat = parameters["covtreat"];
  double curetime = parameters["curetime"];
  
  
  
  //define model compartments
  double S = state["S"];
  double I = state["I"];
  double R = state["R"];
  
  double beta = R0*durinf;
  double omega = 1/durimm;
  

//conditions for theta
  double theta;
  if(t < startvac)
    theta = 0;
  if((startvac <= t) && (t < (startvac+durvac)))
    theta = -(1/durvac)*log(1-(covvac/100));
  if(t > (startvac+durvac))
        theta = 0;
      
//conditions for nu
  double nu;
  if(t < starttreat)
    nu = 1/durinf;
  if(t >= starttreat)
    nu = (1-(covtreat/100))*(1/durinf) + (covtreat/100)*(1/curetime);
          
  double dS = -beta*S*I/P + omega*R - theta*S;
  double dI = beta*S*I/P - nu*I;
  double dR = nu*I - omega*R + theta*S;                    
          
  double prev = 100*I/(S+I+R);
  double pop = S+I+R;
  
  //case per 1000 per year
  double inc = 1000*beta*S*I/P;
  
  //define output parameters
  NumericVector compartments;
  compartments["dS"] = dS;
  compartments["dI"] = dI;
  compartments["dR"] = dR;
  
  List outputlist(4);
  outputlist[0] = compartments; //this variable will be passed into the ODE solver.
  outputlist["prev"] = prev;   //this variable will NOT be passed into the ODE solver.
  outputlist["pop"] = pop;  //this variable will NOT be passed into the ODE solver.
  outputlist["inc"] = inc; //this variable will NOT be passed into the ODE solver.
  return outputlist;
}

