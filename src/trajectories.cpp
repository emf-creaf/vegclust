#include <Rcpp.h>
using namespace Rcpp;

//'
//' Functions for trajectory analysis
//' 
//' Besse, P., Guillouet, B., Loubes, J.-M. & François, R. (2016). Review and perspective for distance based trajectory clustering. IEEE Trans. Intell. Transp. Syst., 17, 3306–3317.
//'
//' De Caceres M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (in preparation). Trajectory analysis in community ecology.
//'
//' Projection of one point into a (directed) segment
//'
//' @param dref Distance between the two segment endpoints
//' @param d1 Distance from the target point to the initial segment endpoint
//' @param d2 Distance from the target point to the final segment endpoint
// [[Rcpp::export(".projectionC")]]
NumericVector projection(double dref, double d1, double d2) {
  double a1 = (pow(d1,2.0)+pow(dref,2.0)-pow(d2,2.0))/(2.0*dref);
  double a2 = dref-a1;
  double s = pow(d1,2.0)-pow(a1,2.0);
  double h = NA_REAL;
  if(s>=0.0) h = sqrt(s);
  return(NumericVector::create(a1,a2,h));
}

//'
//' Distance from one point to one segment
//'
// [[Rcpp::export(".distanceToSegmentC")]]
NumericVector distanceToSegment(double dref, double d1, double d2) {
  NumericVector p = projection(dref,d1, d2);
  if(NumericVector::is_na(p[2]) | (p[0]<0.0) | (p[1]<0.0)) {
    if(d1<d2) {
      p[0] = 0.0;
      p[1] = dref;
      p[2] = d1;
    } else {
      p[0] = dref;
      p[1] = 0.0;
      p[2] = d2;
    }
  }
  return(p);
}

//'
//' Distance between two segments
//'
// [[Rcpp::export(".twoSegmentDistanceC")]]
double twoSegmentDistance(NumericMatrix dmat12, String type="directed-segment") {
  double ds1e1 = dmat12(0,1);
  double ds1s2 = dmat12(0,2);
  double ds1e2 = dmat12(0,3);
  double de1s2 = dmat12(1,2);
  double de1e2 = dmat12(1,3);
  double ds2e2 = dmat12(2,3);
  double Ds = NA_REAL; 
  if((type=="Hausdorff") | (type == "directed-segment")) {
    NumericVector ps1_2 = distanceToSegment(ds2e2,ds1s2, ds1e2);
    NumericVector pe1_2 = distanceToSegment(ds2e2,de1s2, de1e2);
    NumericVector ps2_1 = distanceToSegment(ds1e1,ds1s2, de1s2);
    NumericVector pe2_1 = distanceToSegment(ds1e1,ds1e2, de1e2);
    double ds1_2 = ps1_2[2];
    double de1_2 = pe1_2[2];
    double ds2_1 = ps2_1[2];
    double de2_1 = pe2_1[2];
    //Modifications for directionality of segments
    if(type == "directed-segment") { 
      if(ps1_2[0]>pe1_2[0]) {
        de1_2 = std::min(ds1e1+ps1_2[2], ds1e1+pe1_2[2]);
      }
      if(ps2_1[0]>pe2_1[0]) {
        de2_1 = std::min(ds2e2+ps2_1[2],ds2e2+pe2_1[2]);
      }
    }
    Ds = max(NumericVector::create(ds1_2, de1_2, ds2_1, de2_1));
  } else if (type=="PPA"){ //Perpendicular/Parallel/Angle
    if(ds1e1 > ds2e2) { // Switch roles if longest segment is 1
      ds2e2 = dmat12(0,1);
      ds1s2 = dmat12(0,2);
      de1s2 = dmat12(0,3);
      ds1e2 = dmat12(1,2);
      de1e2 = dmat12(1,3);
      ds1e1 = dmat12(2,3);
    }
    //Assumes longer segment is 2
    NumericVector ps1_2 = distanceToSegment(ds2e2,ds1s2, ds1e2);
    NumericVector pe1_2 = distanceToSegment(ds2e2,de1s2, de1e2);
    double lp1 = ps1_2[2];
    double lp2 = pe1_2[2];
    double dperpendicular = (pow(lp1,2.0)+pow(lp2,2.0))/(lp1+lp2);
    double lpar1 = std::min(ps1_2[0],ps1_2[1]);
    double lpar2 = std::min(pe1_2[0],pe1_2[1]);
    double dparallel = std::min(lpar1,lpar2);
    double dangle = (std::max(lp2,lp1)-std::min(lp2,lp1));
    if(ps1_2[0]>pe1_2[0]) dangle = ds1e1;
          
    Ds = (dperpendicular+dparallel+dangle);
  } 
  return(Ds);
}

//'
//' Triangle inequality for one triplet
//'
// [[Rcpp::export(".triangleinequalityC")]]
bool triangleinequality(double d1, double d2, double d3, double tol=0.0001){
  if((d1+d2)<d3*(1.0-tol)) return(false);
  else if((d1+d3)<d2*(1.0-tol)) return(false);
  else if((d2+d3)<d1*(1.0-tol)) return(false);
  return(true);
}

//'
//' Determines if the matrix fulfills the triangle inequality
//'
// [[Rcpp::export(".ismetricC")]]
bool ismetric(NumericMatrix dmat, double tol=0.0001) {
  int n = dmat.nrow();
  for(int i=0; i<n;i++) {
    for(int j=i;j<n;j++) {
      for(int k=j;k<n; k++) {
        bool ti = triangleinequality(dmat(i,j), dmat(i,k), dmat(j,k), tol);
        if(!ti) return(false);
      }
    }
  }
  return(true);
}

// NOT PRESENTLY USED
double pt(double dIT, double dXT, double dPX, double dPI) {
  if(dPI==0) return(dIT);
  else if(dPX==0) return(dXT);
  double A = pow(dXT,2.0)- pow(dPX,2.0);
  double B = pow(dIT,2.0)- pow(dPI,2.0);
    
  double ax =pow(dXT,2.0)+pow(dIT,2.0);
  double bx = (pow(dXT,2.0)*2.0*B) + (pow(dIT,2.0)*2.0*A) - (4.0*pow(dIT,2.0)*pow(dXT,2.0));
  double cx = pow(dXT,2.0)*pow(B,2.0) + pow(dIT,2.0)*pow(A,2.0);
  double z = pow(bx,2.0)-(4.0*ax*cx);
  double d2 = ((-1.0)*bx + sqrt(z))/(2.0*ax);
  return(sqrt(d2));
}

