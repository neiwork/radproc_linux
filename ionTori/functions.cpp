// Auxiliary parameters
double z1 = 1.0 + pow( 1.0 - (a/mass_bh)*(a/mass_bh), 1.0/3.0) * ( pow(1.0 + a/mass_bh, 1.0/3.0) + 
pow(1.0 - a/mass_bh, 1.0/3.0) );
double z2 = sqrt( 3.0 * (a/mass_bh)*(a/mass_bh) + z1*z1);

// marginally stable circular orbit
double r_ms = mass_bh * (3.0 + z2 - sqrt( (3.0 - z1) * (3.0 + z1 + 2.0*z2) ) );

// marginally bound circular orbit
double r_mb = 2.0 * mass_bh - a + 2.0 * sqrt(mass_bh) * sqrt(mass_bh-a);

// Keplerian specific angular momentum
double angular_mom(double r) {
    return sqrt(mass_bh) * ( r*r - 2.0 * a * sqrt(mass_bh*r) + a*a ) /
               (pow(r, 1.5) - 2.0 * mass_bh * sqrt(r) + a* sqrt(mass_bh) );
}

double l_ms = l_K(r_ms);                                           // Keplerian specific angular momentum at r = r_ms
double l_mb = l_K(r_mb);                                          // Keplerian specific angular momentum at r = r_mb
double l_0  = (1-lambda)*l_ms + lambda*l_mb;          // torus specific angular momentum

 // auxiliary function to find r_cusp and r_center
double aux(double r) { 
    return angular_mom(r) - l_0;
}

//r_cusp = bisect(aux, r_mb, r_ms);      Acá necesito encontrar las dos raíces de aux. Una es el r_cusp
//r_center  = bisect(aux, r_ms, 1000);   y la otra r_center.

/////////////////////////////////////////////////////////////////////////////////////
// METRIC COMPONENTS (in Boyer-Lindquist coordinates)

double g_tt(double r, double theta)  {
    double sigma = r*r + a*a*cos(theta)*cos(theta);
    return -(1.0 - 2.0*mass_bh*r / sigma);
}

double g_rr(double radius, double theta) {
    double delta = r*r - 2.0 * mass_bh * r + a*a;
    double sigma = r*r + a*a*cos(theta)*cos(theta);
    return sigma / delta;
}

double g_thetatheta(double r, double theta) {    // = Sigma
    return r*r + a*a*cos(theta)*cos(theta);
}

double g_tphi(double r, double theta) {
    double sigma = r*r + a*a*cos(theta)*cos(theta);
    return -(2.0*mass_bh*r*a / sigma) * sin(theta)*sin(theta);
}

double g_phiphi(double r, double theta) {
    double sigma = r*r + a*a*cos(theta)*cos(theta);
    return ( r*r + a*a + 2.0*mass_bh*r*a*a*sin(theta)*sin(theta) / sigma ) * sin(theta)*sin(theta);
}

/////////////////////////////////////////////////////////////////////////////////////

// ANGULAR VELOCITIY OF THE TORUS
double angular_vel(double r, double theta)  {
    return - ( g_tphi(r, theta) + l_0 * g_tt(r,theta) ) / ( g_phiphi(r, theta) + l_0 * g_phiphi(r, theta) );
}

// POTENTIAL FUNCTION
double potential(double r, double theta) {
    double aux = g_tt(r,theta) + 2.0*angular_vel(r,theta)*g_tphi(r,theta) + 
    g_phiphi(r,theta) * np.square(angular_vel(r,theta))
    return (aux < 0.0) ? 0.5 * log(-aux / pow(g_tt(r,theta)+angular_vel(r,theta)*g_tphi(r,theta),2 )): 0.0;
}

double potential_s = potential(r_cusp, pi/2.0)         // potential at the torus surface
double potential_c = potential(r_center, pi/2.0)      // potential at the torus center

// NORMALIZED POTENTIAL FUNCTION
double w(double r, double theta) {
    return (potential(r,theta) - potential_s) / (potential_c - potential_s);
}

// ENERGY DENSITY
double energy (double r, double theta) {
    return ( w(r, theta) > 0 ) ? pow(pk, -n) * pow( pow(pk*pow(energy_c, 1./n) 
    + 1.0, w(r,theta)) - 1.0, n) : 0.;
}

// PRESSURE
double pressure(double r, double theta) {
    return pk * pow(energy(r, theta), 1.0 + 1.0/n);
}

/////////////////////////////////////////////////////////////////////////////
// TEMPERATURES

// Electrons
double temp_e(double r, double theta) {
    return (w(r, theta) > 0.0) ? (1.0 - w(r,theta))*M_0 + w(r,theta)*M_1 *
                mu_e * ( (1-beta)*mu*pressure(r,theta) ) / ( kB * energy(r,theta) ) : 0.0;
}

// Ions
double temp_i(double r, double theta) {
    return (w(r, theta) > 0.0) ? ( (1 - w(r,theta))*M_0 + w(r,theta)*M_1 ) *
                mu_i * ( (1-beta)*mu*pressure(r,theta) ) / ( kB * energy(r,theta) ) : 0.0;
}
///////////////////////////////////////////////////////////////////////////