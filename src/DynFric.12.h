#include <math.h>
#include <cmath>
using namespace std;

// Cameron Hummels
// Code for calculating dynamical friction using Chandrasekhar approximation
// Date: 12.2.04
// Version: 12

// This is coded straight from the sigma approximation in Zentner & Bullock
// 2003.  This aproximation is accurate to 1% for x = .01-100.
// In this case, x is the radius of Sgr over the scale radius of the halo
// (I think).  Should only work for calculating sigma for an NFW halo.

double Sigma(double vmax, double x)
{
    double sigma = vmax * 1.4393 * pow(x,0.354) / (1 + 1.1756 * pow(x,0.725));
    return sigma;
}

// This code was written partially through porting KVJ's Fortran code
// and partially through using Zentner & Bullock 2003 as well as Binney
// & Tremaine references, and partially using Hashimoto et al 2003.


void DFAccelNFW(double const msat, double acceleration[3], double const vsat[3],
                double const rsat, double const rhalo, double const GMhalo,
		double const epsilon, double const rsat_scale)

// msat, vsat and rsat are mass, velocity and radius of satellite with respect
// to center of host galaxy (rsat is just the magnitude of the radius, vsat
// is the full vector).  rhalo=d and Gmhalo are from NFW consts.  rhalo
// is scale length for halo.  Gmhalo is total mass of halo times G (grav
// const). acceleration is redefined in this function and passed back
// to the calling function to be added to the gravitational acceleration.

{
    double PI = 3.1415926535;		// Pi
    double G = 0.000004300915; 		// Grav const in sim units
    double G2 = G*G;			// Grav const squared
    double mhalo = GMhalo/G;		// The total mass of the halo
    double vmax = 225.3947;		// maximum velocity of disk in NFW
    double rscale = rsat/rhalo;		// the scaled radius of the sat
    double coulog_hashi, coulog_bullock;

    // the square of the mag of the satellite velocity
    double vsat2 = vsat[0]*vsat[0]+vsat[1]*vsat[1]+vsat[2]*vsat[2];

    // sigma is the one-dimensional velocity dispersion of particles in
    // the host halo.  We use this equation to approximate the Jeans 
    // equation in an NFW profile.  From Zentner & Bullock 2003.
    double sig = Sigma(vmax, rscale);	

    // Chi is ...
    // Defined in both Zentner and KVJ code.
    double chi = sqrt(vsat2)/(sqrt(2.0)*sig);		
    
    // mass of halo interior to radius r -- could be made a function
    // code ported from kvj, from derivation in notes.  specifically
    // for NFW halo
    double mass_int_r = mhalo*(log(rscale+1)-rscale/(rscale+1)); 

    // local density near satellite.  code ported from kvj, from derivation
    // in notes.  specifically for NFW halo
    double rho = mhalo/(4.*PI*rsat*(rsat+rhalo)*(rsat+rhalo));	

    // rtide is tidal radius of satellite? tidal radius of halo?
    // ported from kvj code.  do not understand.
    double rtide = rsat*pow((msat/(3.0*mass_int_r)),(1.0/3.0));

    // coulog is the Coulomb Logarithm
    // determined from Hashimoto 2003 
    // note: coulog becomes 0 if satellite radius is less than 1.4 epsilon
    // this is so we don't get dynamical acceleration at small radii
    if (rsat >= (1.4*epsilon))
    	coulog_hashi = log(rsat / (1.4 * epsilon));
    else
	coulog_hashi = 0.;

    // force_dynfric is the force on the satellite due to dynamical friction.
    // Same equation found in KVJ, Zentner, and Binney & Merrifield.
    // Originally derived from Chandrasekhar.
    double force_dynfric = 4*PI*coulog_hashi*G2*msat*rho*( erf(chi)-2*chi*exp(-(chi*chi))/sqrt(PI))/(vsat2*sqrt(vsat2));


    // the true dynamical friction force opposes the motion of the body
    // by directly opposing the velocity, not the acceleration.
    // Acceleration = force divided by mass of satellite.
    for (int k = 0; k < 3; k++)
	acceleration[k] = -force_dynfric*vsat[k]; 
    return;
}

void DFAccelLog(double const msat, double acceleration[3], double const vsat[3],
    	        double const rsat, double const rhalo, double const v_halo2,
  	        double const epsilon, double const rsat_scale)

// msat, vsat and rsat are mass, velocity and radius of satellite with respect
// to center of host galaxy (rsat is just the magnitude of the radius, vsat
// is the full vector).  rhalo=c and v_halo2 are from log consts.  rhalo
// is scale length for halo.  v_halo2 is square of halo velocity.
// acceleration is redefined in this function and passed back
// to the calling function to be added to the gravitational acceleration.

{
    double PI = 3.1415926535;		// Pi
    double G = 0.000004300915; 		// Grav const in sim units
    double G2 = G*G;			// Grav const squared
    double coulog_hashi, coulog_bullock;

    // the square of the mag of the satellite velocity
    double vsat2 = vsat[0]*vsat[0]+vsat[1]*vsat[1]+vsat[2]*vsat[2];

    // this is the definition of sigma for a logarithmic potential
    double sig = sqrt(v_halo2);

    // Chi is ...
    // Defined in both Zentner and KVJ code.
    double chi = sqrt(vsat2)/(sqrt(2.0)*sig);		
    
    // logarithmic potential-specific mass_int_r
    double mass_int_r = 2 * v_halo2 * rsat / G;

    // logarithmic potential-specific rho
    double rho = 2 * v_halo2 / (4. * PI * G * rsat*rsat);

    // rtide is tidal radius of satellite? tidal radius of halo?
    // ported from kvj code.
    double rtide = rsat*pow((msat/(3.0*mass_int_r)),(1.0/3.0));

    // coulog is the Coulomb Logarithm
    // determined from Hashimoto 2003 
    // note: coulog becomes 0 if satellite radius is less than 1.4 epsilon
    // this is so we don't get dynamical acceleration at small radii
    if (rsat >= (1.4*epsilon))
    	coulog_hashi = log(rsat / (1.4 * epsilon));
    else
	coulog_hashi = 0.;

    // force_dynfric is the force on the satellite due to dynamical friction.
    // Same equation found in KVJ, Zentner, and Binney & Merrifield.
    // Originally derived from Chandrasekhar.
    double force_dynfric = 4*PI*coulog_hashi*G2*msat*rho*( erf(chi)-2*chi*exp(-(chi*chi))/sqrt(PI))/(vsat2*sqrt(vsat2));

    // the true dynamical friction force opposes the motion of the body
    // by directly opposing the velocity, not the acceleration.
    // Acceleration = force divided by mass of satellite.
    for (int k = 0; k < 3; k++)
	acceleration[k] = -force_dynfric*vsat[k]; 
    return;
}
