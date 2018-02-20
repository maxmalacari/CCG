// Max Malacari (KICP - The University of Chicago) - 06/02/2018
// Simulate isotropic or dipolar CRs following the Auger ICRC17 energy spectrum or a single power law
// Boost CRs according to bulk velocity u_motion in direction (l_motion, b_motion)
// Decompose the resulting sky map into spherical harmonics to determine the resulting dipole

// C/C++ classics
#include <iostream>
#include <fstream>
#include <vector>

// Toolkit
#include "maptools.h"
#include "healpixmap.h"
#include "harmotools.h"

// ROOT
#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TVector3.h"

using namespace std;

struct Dipole{
  double amplitude, theta, phi;
};

const string dipoleNames[5] = {"total", "beam", "thresh_mig", "stream", "unboosted"};

void GenerateMonopoleDirections(double & l, double & b);
void GenerateDipoleDirections(double & l, double & b, double theta_dip, double phi_dip, double amp);
double GenThetaFromDipole(double amp);
double GetRandomEnergyICRC(double Emin, double Emax);
double GetRandomEnergyPowerLaw(double min, double index);
double GetRandomEnergyPowerLawBounded(double min, double max, double index);
Dipole CalculateDipole(THealpixMap & inputMap);
void SaveAsFits(THealpixMap & theMap, string fileName);
static inline void loadBar(long long int x, long long int n, int r, int w);

int main(){

  // --- USER SELECTABLE PARAMETERS ---

  const unsigned int nevents = 1000000; // number of random CRs we will simulate

  // Energy related parameters (integral version unless 'differential' is true)
  const bool ICRCspectrum = true; // Use the ICRC17 spectrum? Otherwise power law
  const double spectralIndex = 2.7; // spectral index in power law case
  const double Ethr_low = 1.e19; // eV
  const bool differential = true;
  const double Ethr_hi = 1.1e19; // eV

  // Boost parameters
  const double u_motion = 384; // km/s - only positive (CHECK THIS - I DON'T THINK IT MATTERS), change (l,b) to reverse direction
  const double l_motion = 264.; // deg
  const double b_motion = 48.; // deg

  // Dipole parameters
  const bool genDipole = false; // sample from a dipole? Otherwise a monopole
  const double l_dip = 233.; // galactic longitude of dipole in deg
  const double b_dip = -13.; // galactic latitude of dipole in deg
  const double amp_dip = 0.065; // amplitude of dipole

  // ----------------------------------


  

  const double c = 299792458; // m/s
  const double beta = u_motion * 1.e3 / c;
  const double gamma = 1./sqrt(1. - beta*beta);
  
  if (u_motion*1.e3 > c) {
    cerr << "Error: boost velocity is too high!" << endl;
    return 0;
  }
  
  const double Emin = Ethr_low / (gamma*(1. + fabs(beta))); // min energy to simulate depends on max. possible boost E shift
  double Emax = 1.e21; // eV
  if (differential) Emax = Ethr_hi / (gamma*(1. - fabs(beta)));

  if (differential)
    cout << "[Differential CG calculation]\nEnergy bin = " << Ethr_low << " eV < E < " << Ethr_hi << " eV -> Emin = " << Emin << " eV, Emax = " << Emax << " eV." << endl;
  else if (!differential && !ICRCspectrum)
    cout << "[Integral CG calculation]\nThreshold energy = " << Ethr_low << " eV -> Emin = " << Emin << " eV." << endl;
  else if (!differential && ICRCspectrum)
    cout << "[Integral CG calculation]\nThreshold energy = " << Ethr_low << " eV -> Emin = " << Emin << " eV, (Auger spec.) Emax = " << Emax << " eV." << endl;
  
  TVector3 uVec(1.,1.,1.);
  uVec.SetTheta((90.-b_motion)*TMath::Pi()/180.); // convert to ROOT convention
  uVec.SetPhi((l_motion-180.)*TMath::Pi()/180.);
  uVec.SetMag(1.);

  // Healpix map properties
  const unsigned int nSide = 256;
  const int nPix = 12*nSide*nSide;

  // Solid angle per pixel
  const double solidAnglePerPix = 4.*TMath::Pi() / double(nPix);

  cout << "Simulating bulk motion of u = " << u_motion << " km/s (gamma = " << gamma << ", beta = "
       << beta << "), in direction l = " << l_motion << " deg, b = " << b_motion << " deg." << endl;
  cout << "Generating and boosting events..." << endl;

  // Create output maps
  const int nMaps = 5;
  THealpixMap outputMap[nMaps];
  for (int i=0; i<nMaps; i++) {
    outputMap[i] = THealpixMap(nSide, 'G');
  }

  double l, b, energy, boostedEnergy, angle2Motion, cosTheta;
  int eventCounter = 0, eventCounterECut=0;
  
  while (eventCounter < nevents) {

    if (genDipole) GenerateDipoleDirections(l, b, b_dip, l_dip, amp_dip);
    else GenerateMonopoleDirections(l, b);

    // Now get random energy, perform boost, make energy cut
    if (ICRCspectrum) energy = GetRandomEnergyICRC(Emin, Emax);
    else if (!ICRCspectrum && !differential) energy = GetRandomEnergyPowerLaw(Emin, spectralIndex);
    else if (!ICRCspectrum && differential) energy = GetRandomEnergyPowerLawBounded(Emin, Emax, spectralIndex);
    
    angle2Motion = AngularDistance(l_motion, b_motion, l, b);
    cosTheta = cos(angle2Motion*TMath::Pi()/180.);
    //double boostedEnergyTest = gamma * energy * (1. + beta * cosTheta); // test without decomposition into momentum components

    TVector3 pVec(1.,1.,1.);
    pVec.SetTheta((90.-b)*TMath::Pi()/180.); // convert to ROOT convention
    pVec.SetPhi((l-180.)*TMath::Pi()/180.);
    pVec.SetMag(1.);

    double pDotu = pVec.Dot(uVec); // component of momentum in direction of boost
    double pBoost = (energy/c) * pDotu; // momentum component in direction of boost
    TVector3 pRemainder = (energy/c)*pVec - pBoost*uVec; // remainder of momentum vector perp. to boost direction
    
    double pBoostPrime = gamma * (pDotu + beta) * (energy/c); // boosted momentum component in direction of boost
    TVector3 boostedP = pBoostPrime*uVec + pRemainder; // new momentum vector, boosted component + perp. component
    boostedEnergy = c * boostedP.Mag(); // new energy
    
    double l_boost = boostedP.Phi() * 180./TMath::Pi() + 180.;
    double b_boost = -1. * (boostedP.Theta() * 180./TMath::Pi() - 90.); // convert back to healpix convention

    int ip_boost = outputMap[0].Ip(l_boost, b_boost);
    int ip = outputMap[0].Ip(l, b); // pixel number for original unboosted direction
   
    if (!differential && boostedEnergy >= Ethr_low){
      outputMap[0][ip_boost] += (1. + beta*cosTheta) / solidAnglePerPix; // total map
      outputMap[2][ip] += 1. / solidAnglePerPix; // energy threshold migration map
      eventCounter++;
    }
    if (differential && boostedEnergy >= Ethr_low && boostedEnergy <= Ethr_hi){
      outputMap[0][ip_boost] += (1. + beta*cosTheta) / solidAnglePerPix; // total map
      outputMap[2][ip] += 1. / solidAnglePerPix; // energy threshold migration map
      eventCounter++;
    }
    // Make sure map statistics are in line with the other two maps - don't make energy cut (effect is energy indep.)
    if (eventCounterECut < nevents){
      outputMap[1][ip_boost] += 1. / solidAnglePerPix; // relativistic beaming dipole - energy cut irrelevant
      outputMap[3][ip] += (1. + beta*cosTheta) / solidAnglePerPix; // streaming dipole - energy cut irrelevant
      outputMap[4][ip] += 1. / solidAnglePerPix; // unboosted map
      eventCounterECut++;
    }

    loadBar(eventCounter, nevents, 100, 50);

  }

  
  // Calculate dipole for output maps and save maps
  for (int i=0; i<nMaps; i++){

    // Calculate dipole
    Dipole theDipole = CalculateDipole(outputMap[i]);
    cout << "Dipole [" << dipoleNames[i] << "] is " << theDipole.amplitude << "% in the direction (l,b) = (" << theDipole.phi << ", " << theDipole.theta << ")" << endl;

    // Save maps
    stringstream fileName;
    string mapPrefix;
    if (ICRCspectrum) mapPrefix = "Auger";
    else {
      stringstream isoPref;
      isoPref << "pLaw_" << spectralIndex;
      mapPrefix = isoPref.str();
    }
    string distrib;
    if (genDipole) distrib = "dip";
    else distrib = "iso";
    string type;
    if (differential) type = "diff";
    else type = "int";
    fileName << "./Output/" << mapPrefix << "_" << distrib << "_" << type << "_" << nevents << "_" << dipoleNames[i] << ".fits";
    SaveAsFits(outputMap[i], fileName.str());
    
  }
  
  return 0;
}


void GenerateMonopoleDirections(double & l, double & b){

  static TRandom3 randGen(time(NULL)); // init. random seed on first run
  static bool init = false;
  if (!init){
    cout << "[Sampling directions from monopole]" << endl;
    init = true;
  }
  
  l = randGen.Uniform(0.,360.);
  b = 90. - acos(randGen.Uniform(-1,1))*180./TMath::Pi();
  
}


void GenerateDipoleDirections(double & l, double & b, double theta_dip, double phi_dip, double amp){

  static TRandom3 randGen(time(NULL)); // init. random seed on first run
  static bool init = false;
  static TVector3 dipoleDir(1.,1.,1.), dir(1.,1.,1.);
  static double theta_dip_corr = TMath::Pi()/180. * (90. - theta_dip);
  static double phi_dip_corr = phi_dip*TMath::Pi()/180.;
  if (!init){
    cout << "[Sampling directions from dipole]: (l, b) = (" << phi_dip << ", " << theta_dip << "), d = " << amp*100. << "%." << endl;
    dipoleDir.SetTheta(theta_dip_corr);
    dipoleDir.SetPhi(phi_dip_corr);
    dipoleDir.SetMag(1.);
    dir.SetMag(1.);
    init = true;
  }

  // Get random directions with dipole oriented along z-axis
  double phi0 = randGen.Uniform(0.,2.*TMath::Pi());
  double theta0 = GenThetaFromDipole(amp);

  dir.SetTheta(theta0);
  dir.SetPhi(phi0);

  // Rotate to true dipole direction
  dir.RotateUz(dipoleDir);
  
  l = dir.Phi()*180./TMath::Pi();
  b = 90. - dir.Theta()*180./TMath::Pi();
  
}

double GenThetaFromDipole(double amp){

  static TRandom3 randGen(time(NULL));

  bool accept = false;

  double theta;
  while (!accept){
    theta = acos(-1.+ 2.*randGen.Uniform(0,1));
    double x = (1.+amp)*randGen.Uniform(0,1);
    if (x <= 1.+amp*cos(theta)) accept = true;
  }

  return theta;
    
}


double GetRandomEnergyICRC(double Emin, double Emax){

  static TF1 spectrum;
  static bool init = false;

  if (!init) {

    cout << "[Sampling energies from Auger ICRC17 energy spectrum]" << endl;
    
    // Combined spectrum fit parameters ICRC17
    const double J0 = 2.78e-19; 
    const double E_ankle = 5.08e18;
    const double E_s = 3.9e19;
    const double gamma1 = 3.293;
    const double gamma2 = 2.53;
    const double deltaGamma = 2.5;
  
    stringstream funcStr1, funcStr2;
    funcStr1 << "(x<=" << E_ankle << ")*" << J0 << "*(x/" << E_ankle << ")^" << -1.*gamma1;
    funcStr2 << "+(x>" << E_ankle << ")*" << J0 << "*(x/" << E_ankle << ")^" << -1.*gamma2 << "*(1+(" << E_ankle/E_s << ")^" << deltaGamma << ")*(1+(x/" << E_s << ")^" << deltaGamma << ")^-1";

    string specFunc = funcStr1.str() + funcStr2.str();
    
    //spectrum = TF1("spectrum",specFunc.c_str(),3.e17,1.e20); // for plotting only. DO NOT USE FOR SIMS
    spectrum = TF1("spectrum",specFunc.c_str(),Emin,Emax);
    spectrum.SetNpx(1000000);

    // TCanvas specCanvas;
    // specCanvas.SetLogx();
    // specCanvas.SetLogy();
    // spectrum.GetXaxis()->SetTitle("Energy [eV]");
    // spectrum.GetYaxis()->SetTitle("J(E) [eV^{-1} km^{-2} sr^{-1} yr^{-1}]");
    // spectrum.GetXaxis()->CenterTitle();
    // spectrum.GetYaxis()->CenterTitle();
    // spectrum.SetLineColor(kAzure+7);
    // spectrum.Draw();
    // specCanvas.SaveAs("./Output/spectrum.pdf");
    
    init = true;
  }

  
  return spectrum.GetRandom();
}

// Sample a power law to infinity
double GetRandomEnergyPowerLaw(double min, double index){

  static bool init = false;
  if (!init){
    cout << "[Sampling energies from power law with gamma = " << index << "]" << endl;
    init = true;
  }

  static TRandom3 randGen(time(NULL));
  double dice = randGen.Uniform(0.,1.);

  return min * pow(1.-dice,1./(-index+1));
}

// Simple single power law
double GetRandomEnergyPowerLawBounded(double Emin, double Emax, double index){

  static TF1 spectrum;
  static bool init = false;
  
  if  (!init) {
    cout << "[Sampling energies from bounded power law with gamma = " << index << "]" << endl;
    stringstream funcStr;
    funcStr << "x^(-1.*" << index << ")";
    
    spectrum = TF1("spectrum",funcStr.str().c_str(),Emin,Emax);

    spectrum.SetNpx(1000000);
    
    init = true;
  }
  
  return spectrum.GetRandom();
}


Dipole CalculateDipole(THealpixMap & inputMap){

  static const double lmax = 1;

  Dipole theDipole;

  // Calculate alms
  
  vector<vector<complex<double> > > alms =  anafast_alm_healpix(inputMap, lmax);
  double a00 = alms[0][0].real();
  double a10 = alms[1][0].real();
  double a11_re = alms[1][1].real();
  double a11_im = alms[1][1].imag();

  // dipole components
  double dx = sqrt(4.*TMath::Pi())*(-sqrt(3./2./TMath::Pi())*a11_re/a00);
  double dy = sqrt(4.*TMath::Pi())*( sqrt(3./2./TMath::Pi())*a11_im/a00);
  double dz = sqrt(4.*TMath::Pi())*(0.5*sqrt(3./TMath::Pi())*a10/a00);

  // dipole amplitude
  double dipole_amp = sqrt(dx*dx + dy*dy + dz*dz);

  // dipole direction
  double dipole_theta = acos(dz/dipole_amp);
  double dipole_phi   = atan2(dy, dx);

  theDipole.amplitude = dipole_amp*100.;
  theDipole.theta = -1. * (dipole_theta*180./TMath::Pi() - 90.); // Coordinates consistent with input
  theDipole.phi = fmod(dipole_phi*180./TMath::Pi() + 360.,360.);

  return theDipole;
}


// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
static inline void loadBar(long long int x, long long int n, int r, int w)
{
    // Only update r times.
    if ( x % (n/r +1) != 0 ) return;
 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++)
       printf("=");
 
    for (int x=c; x<w; x++)
       printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n\033[F\033[J");
}

// Silly function to save healpix map as a fits file
void SaveAsFits(THealpixMap & theMap, string fileName){
  char* chrName = const_cast<char*>(fileName.c_str());
  theMap.WriteFits(chrName);
}
