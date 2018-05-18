//saralamba@gmail.com
//This is an example of using Rcpp with deSolve.
// After you compile this file with the command Rcpp::sourceCpp, 
// the function SIRmodel will be available in the Global environment.


#include <Rcpp.h>
//#include <math.h>  uncomment this line if you have the math functions, i.e. cos, tan, and log in your model
using namespace Rcpp;


// [[Rcpp::export]]
List SIRmodel(double t, NumericVector state, NumericVector parameters) {

  //define parameters
  double birth = parameters["birth"];
  double beta = parameters["beta"];
  double death = parameters["death"];
  double recovery = parameters["recovery"];
  
  //define model compartments
  double S = state["S"];
  double I = state["I"];
  double R = state["R"];
  
  //define model structure
  double dS = birth - beta*I*S - death*S;
  double dI = beta*I*S - recovery*I - death*I;
  double dR = recovery*I - death*R;
  
  double extraParm = beta*I*S;
  
  //define output parameters
  NumericVector compartments;
  compartments["dS"] = dS;
  compartments["dI"] = dI;
  compartments["dR"] = dR;
  
  List outputlist(2);
  outputlist[0] = compartments; //this variable will be passed into the ODE solver.
  outputlist["extraParm"] = extraParm;   //this variable will NOT be passed into the ODE solver.
  
  return outputlist;
}

