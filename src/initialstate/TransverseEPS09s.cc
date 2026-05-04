#include "TransverseEPS09s.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>

//Constructor
fit_parameters::fit_parameters(int eps_order, int eps_pset){
  order = eps_order, pset = eps_pset;
  for(int n = 0;n < 31;n++){
    std::string filename = get_filename(order,n+1);
    std::ifstream infile(filename.c_str());
    if(!(infile.good())){
      std::cout << "No file " << filename << " found!" << std::endl;
      std::cout << "Exiting program..." << std::endl;
      infile.close();
      exit(1);
    }
    std::string dummy;
    std::string file_version;
    infile >> dummy;
    infile >> file_version;
    if(file_version != "1.1"){ // Checks that file is valid
      std::cout << "File: " << filename << " not valid!" << std::endl;
      std::cout << "Exiting program..." << std::endl;
      exit(1);
    }

    //Read parameter values from file:
    for(unsigned int i = 0;i < 8;i++){ //Flavors    
      for(unsigned int j = 0;j < 51;j++){ //Q-grid
	for(unsigned int k = 0;k < 51;k++){ //x-grid
	  for(unsigned int l = 0; l < 4;l++){ //c_i
	    infile >> c_values[n][i][j][k][l];
	  }
	}
      }
    }
    infile.close();
  }
}

/*
 * Returns the static instance of fit parameters (singleton pattern).
 * Reads new files if different order is requested.
 */
fit_parameters* fit_parameters::_instance = 0;
fit_parameters* fit_parameters::getInstance(int eps_order, int eps_pset){
  if(_instance == 0)
    _instance = new fit_parameters(eps_order, eps_pset);
  if(_instance->order == eps_order){
    return _instance;    
  } else {
    _instance = new fit_parameters(eps_order, eps_pset);
    return _instance;
  }
}

/*
 * Constructs the file name for selected order & parameter set
 * using string stream from STD library
 *
 * "eps09s" + ORDER + PSET + ".dat"
 * ORDER = LO, NLO  PSET = 1,...,31
 */
std::string fit_parameters::get_filename(int eps_order, int eps_pset){
  std::string filename;
  switch(eps_order){
  case 1:
    filename = "../examples/eps09s/eps09sLO";
    break;
  case 2:
    filename = "../examples/eps09s/eps09sNLO";
    break;
  default: //Checks that order is valid
    std::cout << "Invalid EPS09s order!! order = " << eps_order << std::endl;
    std::cout << "Should be 1 or 2" << std::endl;
    std::cout << "Exiting program..." << std::endl;
    filename = " ";
    exit(2);
    break;
  }  
  if( (eps_pset < 1) || (eps_pset > 31) ){ //Checks that pset is valid
    std::cout << "Invalid EPS09s parameters set!! set = " << eps_pset << std::endl;
    std::cout << "Should be between 1 and 31" << std::endl;
    std::cout << "Exiting program..." << std::endl;
    exit(3);
  }
  std::stringstream ss_pset;
  ss_pset << eps_pset;
  filename += ss_pset.str();
  filename += ".dat";
  return filename;
}

/*
 * Polynomial interpolation with Newton's divided difference method.
 */
double pol_int(double f_i[], double x_i[], int n, double x){
  for(int i = 1;i < n;i++){
    for(int j = n-1;j > i - 1;j--){
      f_i[j] = (f_i[j] - f_i[j-1])/(x_i[j] - x_i[j-i]);
    }
  }

  double p = f_i[n-1];
  for(int i = n-2;i > -1;i--){
    p = (x - x_i[i])*p + f_i[i];
  }
  
  return p;
}

/*
 * Interpolation for fit parameters for given x, q for all flavors
 */
void fit_parameters::interpolate_c_values(const double xx, const double q, 
					  const int pset, double c[8][4]){
  double xmin = 1e-6, xmax = 1, q2min = 1.69, q2max = 1000000;
  double xcut = 0.1;
  int nx = 51, nq = 51; 
  int nxlog = 25, nxlin = 25;
  int nparams = 4;

  //Freeze Q^2 if outside the limits
  double q2 = q*q;
  if(q2 > q2max){
    q2 = q2max;
  } else if (q2 < q2min){
    q2 = q2min;
  }

  //Freeze x values if outside the limits
  double x = xx;
  if(x > xmax){
    x = xmax;
  } else if (x < xmin){
    x = xmin;
  }

  //Calculate the position in log(log Q^2) grid:
  double i_qd = (nq-1)*log( log(q2)/log(q2min) )/log( log(q2max)/log(q2min) );
  int i_q = static_cast<int>( i_qd );

  //Set the q-index to interval [1,...,49]
  if(i_q < 1){
    i_q = 1;
  } else if (i_q > (nq-2) ){
    i_q = nq-2;
  }

  //Calculate the three nearest points in log(log Q^2) grid
  double q_i[3];
  for(int i = 0;i < 3;i++){
    q_i[i] = sqrt( exp( pow(log(q2max),(i_q+i-1)*1.0/(nq-1))*
			pow(log(q2min),1-(i_q+i-1)*1.0/(nq-1)) ) );
  }

  //Calculate the position in log(x) or x grid
  int i_x;
  if(x <= xcut){
    i_x = static_cast<int>( nxlog*log(x/xmin)/log(xcut/xmin) );
  } else {
    i_x = static_cast<int>( (x-xcut)*nxlin/(xmax-xcut) + nxlog );
  }

  //Set the x-index to interval [1,...,48]
  if(i_x < 1){
    i_x = 1;
  } else if (i_x > nx - 3){
    i_x = nx - 3;
  }

  //Calculate the three nearest points in log(x) or x grid
  double x_i[4];
  for(int i = 0;i < 4;i++){
    if(i_x-1+i < nxlog){
      x_i[i] = xmin*exp( ((i_x-1+i)*1.0/nxlog)*log(xcut/xmin) );
    } else {
      x_i[i] = ( ( i_x - 1 + i - nxlog)*1.0/nxlin )*(xmax-xcut) + xcut;
    }
  }

  //Interpolate c parameters:
  double ccc[8][nparams][3][4], cc[8][nparams][3];
  for(int i = 0; i < 8;i++){
    for(int j = 0; j < nparams;j++){
      for(int k = 0; k < 3;k++){
	for(int l = 0; l < 4;l++){
	  ccc[i][j][k][l] = c_values[pset-1][i][i_q+k-1][i_x+l-1][j];
	}
	cc[i][j][k] = pol_int(ccc[i][j][k], x_i, 4, x);
      }
      c[i][j] = pol_int(cc[i][j], q_i, 3, sqrt(q2));
    }
  }
  return;
}


/*
 * Recursive function for adaptiveBoole function. Splits the interval so that
 * the required accuracy is obtained
 */ 
double adaptiveBooleRec(double (*f)(double, void*), void *p,
			double x1, double x9, double epsilon, double S,
			double f1, double f3, double f5, double f7, 
			double f9, int rec_steps){
  double h = x9 - x1, x5 = x1 + h/2;
  double x2 = x1 + h/8, x4 = x1 + 3*h/8, x6 = x1 + 5*h/8, x8 = x1 + 7*h/8;
  double f2 = f(x2,p), f4 = f(x4,p), f6 = f(x6,p), f8 = f(x8,p);
  double S_left = (h/180)*(7*f1 + 32*f2 + 12*f3 + 32*f4 + 7*f5);
  double S_right = (h/180)*(7*f5 + 32*f6 + 12*f7 + 32*f8 + 7*f9);
  double S_tot = S_left + S_right;
  if (rec_steps <= 0 || fabs(S_tot - S) <= 15*epsilon){
    return S_tot + (S_tot - S)/15;
  }
  return adaptiveBooleRec(f,p,x1,x5,epsilon/2, S_left, f1, f2, f3, f4, f5,
			  rec_steps-1)
    +    adaptiveBooleRec(f,p,x5,x9,epsilon/2, S_right, f5, f6, f7, f8, f9,
			  rec_steps-1);
}         
 
/*
 * Adaptive Boole's Rule
 *
 * Input:
 *  double (*f)(double,void*) = pointer to function with unknown parameters
 *  double x1, x5             = interval [a,b]
 *  double epsilon            = error tolerance
 *  int max_recursion_steps     = maximal number of recursions
 *
 * Output:
 *  the result from the integration (double)
 *
 * Parameters tested to be suitable for this purpose. 
 * For other purposes use with caution (or not at all)!
 */ 
double adaptiveBoole(double (*f)(double, void* ), void *p, 
	       double x1, double x5, double epsilon, int max_recursion_steps){
  double h = x5 - x1, x2 = x1 + h/4, x3 = x1 + h/2, x4 = x1 + 3*h/4;
  double f1 = f(x1,p), f2 = f(x2,p), f3 = f(x3,p), f4 = f(x4,p), f5 = f(x5,p);
  double S = (h/90)*(7*f1 + 32*f2 + 12*f3 + 32*f4 + 7*f5);
  return adaptiveBooleRec(f,p,x1,x5,epsilon, S, f1, f2, f3, f4, f5,
			  max_recursion_steps);
}                   

struct paramsi1d1 { int i1; double d1; };

/*
 * Calculates the Woods-Saxon distribution for given position.
 *
 * Mapped with change of variables as z = (1-t)/t so that the integral from
 * 0 < z < Infinity becomes an integral where 0 < t < 1.
 * Woods-Saxon parameters are 
 *
 *  R_A = 1.12*A^(1/3) - 0.86*A(-1/3) fm
 *    d = 0.54 fm
 *  n_0 = 3*A/(4*Pi*R_A^3)*1/(1 + (Pi*d/R_A)^2) (Correct normalization for A>3)
 */
double mapped_woodSaxon_density(double t, void *params){
  if(t == 0) return 0;
  paramsi1d1 * wsparams = (paramsi1d1*)params;
  int a = wsparams->i1;
  double s = wsparams->d1;
  double d = 0.54;
  double r = 1.12*pow(a,1./3.) - 0.86*pow(a,-1./3.);
  double n0 = 0.75*a/(acos(-1.0)*pow(r,3))/(1+(pow(acos(-1.0)*d,2))/(r*r));
  return 2*n0/( t*t*(1 + exp( (sqrt( s*s + (1-t)*(1-t)/(t*t) ) - r )/d ) ) );
}

/*
 * Calculates T_A^{WS}(s) in given point s for given A
 */
double taWoodsSaxon(int a, double s){
  paramsi1d1 params = { a, s };
  return adaptiveBoole(mapped_woodSaxon_density, &params, 0,1,1e-12, 12);
}

  /* Deuterium thickness function */

//S-wave 
double deuterium_u(double r, double beta, double gamma, double epsilon, 
		   double xc, double alpha, double n){
  if(alpha*r > xc){
    return n*sqrt(1 - epsilon*epsilon)*( 1 - exp(-beta*(alpha*r-xc) ) )
      *exp(-alpha*r)/r;
  } else {
    return 0;
  }
}
//D-wave
double deuterium_w(double r, double beta, double gamma, double epsilon, 
		  double xc, double alpha, double n){
  if(alpha*r > xc){
    return n*epsilon*( 1 - exp( -gamma*(alpha*r-xc) ) )*
      ( 1 - exp( -gamma*(alpha*r-xc) ) )*exp(-alpha*r)*
      ( 1 + 3*(( 1 - exp(-gamma*alpha*r) )/(alpha*r))*
	( 1 + ( 1 - exp(-gamma*alpha*r) )/(alpha*r) ) )/r;
  } else {
    return 0;
  }
}

/*
 * Calculates the square of the wave function (two possible parameter set from
 * Nucl.Phys. A730 (2004) 448-459)
 */
double deuterium_psi2_ds(double r, int d_set){
  double alpha = 1.0/4.316;
  double beta, gamma, epsilon, xc, rho, n;
  switch(d_set){
  case 1:
    beta = 4.680, gamma = 2.494, epsilon = 0.03232, xc = 0;
    rho = -27.944041219;
    n = sqrt( 2*alpha/(1 - alpha*rho) );
    break;
  case 2:
    beta = 9.045, gamma = 4.799, epsilon = 0.02438, xc = 0.13;
    rho = -28.1136;
    n = sqrt( 2*alpha/(1 - alpha*rho) );
    break;
  default:
    std::cout << "Wrong deuteron wave function parameter set! " << std::endl;
    exit(4);
    break;
  }
  double w = deuterium_w(r, beta, gamma, epsilon, xc, alpha, n);
  double u = deuterium_u(r, beta, gamma, epsilon, xc, alpha, n);
  return (u*u + w*w);
}

/*
 * Converts the wave function from p-n distance to distance from the CM
 */
double deuterium_psi2_ds_cm(double r, int d_set){
  return 8*deuterium_psi2_ds(2*r, d_set);
}

/*
 * Integrand for t_deuterium
 */
double mapped_t_d_integrand(double rl, void *data){
  if(rl == 0) return 0; 
  paramsi1d1 * parameters = (paramsi1d1*)data;
  int paramset = parameters->i1;
  double rt = parameters->d1;
  return 2*deuterium_psi2_ds_cm(sqrt(rt*rt + (1-rl)*(1-rl)/(rl*rl)), 
				paramset)/(rl*rl);
}

/*
 * Deuterium thickness function from Hulthen wave function
 */
double t_deuterium(double rT){
  paramsi1d1 params = { 1, rT };
  return 2*adaptiveBoole(mapped_t_d_integrand, &params, 0,1,1e-12, 12);
}


/*
 * Constructor:
 *
 * Calculates thickness function values for given A for several s values
 */
thickness_function::thickness_function(int aa){
  a = aa;
  n_points_lins = 150;
  n_points_linu = 50;
  tail_length = 4; //[fm]
  double r = 1.12*pow(a,1./3.) - 0.86*pow(a,-1./3.);
  double s_cut = r + tail_length;
  
  if(a == 2){ //Deuterium wave function
    for(int i = 0; i < n_points_lins; i++){
      double s = i*(s_cut)/(1.0*n_points_lins);
      ta[i] = t_deuterium(s);
    }
    for(int i = 0; i < n_points_linu; i++){
      double s = s_cut + s_cut/n_points_lins*i*i;
      ta[i+n_points_lins] = t_deuterium(s);
    }
  } else if (a > 2){ //Woods-Saxon
    for(int i = 0; i < n_points_lins; i++){
      double s = i*(s_cut)/(1.0*n_points_lins);
      ta[i] = taWoodsSaxon( a, s );
    }    
    for(int i = 0; i < n_points_linu; i++){
      double s = s_cut + s_cut/n_points_lins*i*i;
      ta[i+n_points_lins] = taWoodsSaxon( a, s );
    }
  } else { //No thickness for A < 2, Program stops!
    std::cout << "Error: no thickness function for A < 2!" << std::endl;
    exit(5);
  }
}

/*
 * Thickness function interpolation:
 *  Interpolates with linear interval in s from 0 to s = R_A + 4fm and 
 *  from this on with linear interval in u = 1/(1 + s).
 */
double thickness_function::operator()(double s){  
  s = fabs(s);
  //std::cout<<" Calculating thickness function for s = " << s << " fm" << std::endl;
  int ns = this->n_points_lins;
  int nu = this->n_points_linu;
  double aa = this->a;
  double tl = this->tail_length;
  double r = 1.12*pow(aa,1./3.) - 0.86*pow(aa,-1./3.);
  double s_cut = r + tl;
  int index1, index2;
  double s1, s2, interpolated_ta;
  if( s < r + tl ){
    index1 = std::min(static_cast<int>(floor( ns*s/s_cut )), ns - 1);
    index2 = index1 + 1;

    s1 = index1*s_cut/ns;
    s2 = index2*s_cut/ns;

    if(index1 != index2){
      interpolated_ta = ta[index1]*exp( log(ta[index2]/ta[index1])
      				*(s-s1)/(s2-s1) );
    } else {
      interpolated_ta = ta[index1];
    }

  } else {
    index1 = std::min(static_cast<int>(floor( sqrt( (ns/s_cut)*(s - s_cut) ) )),
		 nu - 1);
    index2 = std::min(index1 + 1, nu - 1 );

    s1 = s_cut + s_cut*index1*index1/ns;
    s2 = s_cut + s_cut*index2*index2/ns;

    if(index1 == index2){
      interpolated_ta = 0;
    } else {
      interpolated_ta = ta[index1 + ns]*exp( log(ta[index2+ns]/ta[index1+ns])
					     *(s-s1)/(s2-s1) );
    }
  }  

  return interpolated_ta;
}

/*
 * Return an instance of static thickness function object (singleton pattern)
 * Can handle thickness functions for two different nucleus.
 */
thickness_function* thickness_function::_instancea = 0;
thickness_function* thickness_function::_instanceb = 0;
thickness_function* thickness_function::getInstance(int aa){
  if(_instancea == 0){
    _instancea = new thickness_function(aa);
  }
  if(_instancea->a == aa){
    return _instancea;
  } else {
    if(_instanceb == 0){
      _instanceb = new thickness_function(aa);
    }
    if(_instanceb->a == aa){
      return _instanceb;
    } else {
      std::cout << "Support only for 2 nucleus in a run! Exiting program..." << std::endl;
      exit(6);
    }
  }
}

 //End of eps09_s namespace


// Initialize EPS09 nPDFs with given order (1=LO, 2=NLO) and error set.

void EPS09s::init(int iOrderIn, int iSetIn) {
  // Save the order and error set number.
  iOrder = iOrderIn;
  iSet   = iSetIn;

  //std::cout<<"s init "<<PDF::sNowA<<" "<<PDF::sNowB<<std::endl;
  double a = getA();
  //std::cout << "Initializing EPS09s for A = " << a << std::endl;
  if(a > 2){
  if(a < 16){
      std::cout << "EPS09s Warning:" << std::endl;
      std::cout << "Only nuclei with A >= 16 used for fitting!"
	   << std::endl;
  }

  ta = thickness_function::getInstance(a);
  fit_params = fit_parameters::getInstance(iOrder, iSet);

}
}

//--------------------------------------------------------------------------

// Interpolation from the grid.

void EPS09s::rUpdate(int , double x, double Q2) {
  //std::cout << "rUpdate this = " << this << std::endl;
  double tas = (*ta)(Tpos);
  //std::cout << "PDF side = " << sideLabel
  //        << " this = " << this
  //        << " Tpos = " << Tpos << std::endl;
  double c[8][4];
  fit_params->interpolate_c_values( x, sqrt(Q2), iSet, c);
  double r[8];
  for(int i = 0;i < 8;i++){
    //r[i] = tas;
    r[i] = 1;
    for (int j = 0;j < 4;j++){
	    //r[i] += c[i][j]*pow(tas,j+2);
      r[i] += c[i][j]*pow(tas,j+1);
    }
  }
  ruv = r[0], rdv = r[1], ru = r[2], rd = r[3], rs = r[4], rc = r[5], 
    rb = r[6], rg = r[7];
}