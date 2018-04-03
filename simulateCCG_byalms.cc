// Max Malacari (KICP - The University of Chicago) - 02/04/2018
// Create a skymap using alms and boost that skymap

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
#include "TH1.h"
#include "TCanvas.h"

using namespace std;

// Structure to hold dipole information
struct Dipole{
  double amplitude, theta, phi;
};

Dipole CalculateDipole(THealpixMap & inputMap);
void SaveAsFits(THealpixMap & theMap, string fileName);

int main(){

  // --- USER SELECTABLE PARAMETERS ---

  const double spectralIndex = 2.7; // spectral index in power law case
  // Boost parameters
  const double u_motion = 384.; // km/s - only positive (CHECK THIS - I DON'T THINK IT MATTERS), change (l,b) to reverse direction
  const double l_motion = 264.; // deg
  const double b_motion = 48.; // deg

  // ----------------------------------

  const double c = 299792458; // m/s
  const double beta = u_motion * 1.e3 / c;
  
  TVector3 uVec(1.,1.,1.);
  uVec.SetTheta((90.-b_motion)*TMath::Pi()/180.); // convert to ROOT convention
  uVec.SetPhi((l_motion-180.)*TMath::Pi()/180.);
  uVec.SetMag(1.);

  // Healpix map properties
  const unsigned int nSide = 256;
  const int nPix = 12*nSide*nSide;

  // Solid angle per pixel
  const double solidAnglePerPix = 4.*TMath::Pi() / double(nPix);

  // Choose some alms to make the map
  double a00 = 1.;
  double a10 = 0.;
  double a20 = 0.;
  double a11_re = 0.;
  double a11_im = 0.;
  double a21_re = 0.;
  double a21_im = 0.;
  double a22_re = 0.;
  double a22_im = 0.;
  vector<vector<complex<double> > > alms_mon = { {{a00,0}} };  
  vector<vector<complex<double> > > alms_dip = { {{a00,0}}, {{a10,0},{a11_re,a11_im}} };
  vector<vector<complex<double> > > alms_quad = { {{a00,0}}, {{0,0},{0,0}}, {{a20,0},{a21_re,a21_im},{a22_re,a22_im}} };
  vector<vector<complex<double> > > alms_comb = { {{a00,0}}, {{a10,0},{a11_re,a11_im}}, {{a20,0},{a21_re,a21_im},{a22_re,a22_im}} };

  // Form the map
  THealpixMap inMap(nSide, 'G');
  inMap.Alm2Map(alms_comb);

  // Try distorting the original sky map using the pixel totals (rather than individual events)
  THealpixMap boostedMap(nSide, 'G');
  // Loop over pixels in the map
  for (int pix=0; pix<nPix; pix++){
    double l,b,val;
    inMap.GiveLB(pix, l, b);
    //val = outputMap[4][pix];
    val = inMap[pix];
    
    // Calculate val shift
    double ang2Boost = AngularDistance(l_motion, b_motion, l, b);
    double cosTheta = cos(ang2Boost*TMath::Pi()/180.);
    double valNew = val*(1. + (spectralIndex - 1.)*beta*cosTheta) * (1. + beta*cosTheta) * (1. + beta*cosTheta) * (1. + beta*cosTheta); // energy shift and streaming

    // Create a TVector3 corresponding to this pixel
    TVector3 pixVec(1.,1.,1.);
    pixVec.SetTheta((90.-b)*TMath::Pi()/180.); // convert to ROOT convention
    pixVec.SetPhi((l-180.)*TMath::Pi()/180.);
    pixVec.SetMag(1.);

    TVector3 rotationAxis = pixVec.Cross(uVec); // Rotate the vector about the axis perpendicular to the boost direction and the pixel
    double rotationAngle = acos(cosTheta) - acos((cosTheta+beta)/(1.+beta*cosTheta)); // Rotation angle due to rel. aberration
    //cout << rotationAngle << endl;
    pixVec.Rotate(rotationAngle, rotationAxis);

    // Get new l and b
    double l_new = pixVec.Phi() * 180./TMath::Pi() + 180.;
    double b_new = -1. * (pixVec.Theta() * 180./TMath::Pi() - 90.); // convert back to healpix convention

    // Fill the map
    int ip = boostedMap.Ip(l_new, b_new);
    boostedMap[ip] += valNew;

    //cout << l << " " << b << " " << val << " " << l_new << " " << b_new << " " << valNew << endl;
  }
  string outname = "./Output/outMap_byalms.fits";
  SaveAsFits(boostedMap, outname);
  Dipole theDip = CalculateDipole(boostedMap);
  cout << theDip.amplitude << endl;
  
  return 0;
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

// Silly function to save healpix map as a fits file
void SaveAsFits(THealpixMap & theMap, string fileName){
  char* chrName = const_cast<char*>(fileName.c_str());
  theMap.WriteFits(chrName);
}


