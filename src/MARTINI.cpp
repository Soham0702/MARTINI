// MARTINI.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009-2010 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the main class
// MARTINI: provide the main user interface to everything else.
// it also contains the main evolution routine

#define AMYpCut 4.01
#include "LorentzBooster.cxx"
#include "MARTINI.h"
#include <iomanip>
#include "LHAPDF/LHAPDF.h"

// Average-nucleon nuclear PDF: (Z/A)*f^{p/A} + (N/A)*f^{n/A},
// with f^{n/A} from f^{p/A} by isospin (u<->d). Reproduces the ensemble
// average of the old per-event next(n1,n2) p/n sampling, deterministically.
class NuclearAveragePDF : public Pythia8::PDF {
public:
  NuclearAveragePDF(std::string setName, int member, int Zin, int Ain)
    : Pythia8::PDF(2212), Z(Zin), A(Ain) {
    lha = LHAPDF::mkPDF(setName, member);
  }
  ~NuclearAveragePDF() override { delete lha; }
private:
  LHAPDF::PDF* lha;
  int Z, A;
  void xfUpdate(int, double x, double Q2) override {
    double xx = std::max(x,  lha->xMin());
    double q2 = std::min(std::max(Q2, lha->q2Min()), lha->q2Max());
    double zf = double(Z)/double(A);   // proton fraction
    double nf = 1.0 - zf;              // neutron fraction
    // bound-proton flavors from the EPPS21 grid (returns x*f)
    double up  = lha->xfxQ2( 2, xx, q2), dp  = lha->xfxQ2( 1, xx, q2);
    double ubp = lha->xfxQ2(-2, xx, q2), dbp = lha->xfxQ2(-1, xx, q2);
    // isospin: u_n = d_p, d_n = u_p (same for sea)
    xu    = zf*up  + nf*dp;
    xd    = zf*dp  + nf*up;
    xubar = zf*ubp + nf*dbp;
    xdbar = zf*dbp + nf*ubp;
    xs    = lha->xfxQ2( 3, xx, q2);
    xsbar = lha->xfxQ2(-3, xx, q2);
    xc    = lha->xfxQ2( 4, xx, q2);
    xb    = lha->xfxQ2( 5, xx, q2);
    xg    = lha->xfxQ2(21, xx, q2);
    xgamma = 0.;
    xcbar = xc;
    xbbar = xb;
  }
};




// constructor
MARTINI::MARTINI()
{
    import  = new Import();
    importLRates  = new ImportLRates();
    random  = new Random();
    testing = new Testing();
    hydroSetup = new HydroSetup();
    hqRates = new HQRates();
    glauber = new Glauber();
    settings.init();
    gsl_rand = gsl_rng_alloc(gsl_rng_taus);
    binary_info_ptr = new binarycollision_info();
}

// destructor
MARTINI::~MARTINI()
{
    delete import;
    delete importLRates;
    delete random;
    delete rates;
    delete testing;
    delete hydroSetup;
    delete hqRates;
    delete glauber;
    delete gsl_rand;
    delete binary_info_ptr;
}

// evolve every parton by one time step. This is the core of MARTINI.
int MARTINI::evolve(vector<Parton>  *plist, vector<Source> *slist, int counter, int it)
{
    HydroInfo hydroInfo;
    ReturnValue f;                  // will contain \Delta p, the change in momentum due to a process
    Norms n;                        // holds the integrals \int dGamma_i dp, currently three for all radiative process
    double qqRate, qgRate;          // hold the integrals \int dGamma_i domega for the elastic processes
    double gqRate, ggRate;
    double conversionqg;            // hold the total conversion rate q->g
    double conversionqgamma;        // hold the total conversion rate q->gamma
    double conversiongq;            // hold the total conversion rate g->q
    
    Vec4 vecp, newvecp;             // momentum four-vectors
    Vec4 vecpRest, newvecpRest;     // momentum four-vectors in rest frame of fluid cell
    Vec4 vecq;                      // four-vectors of momentum transfer
    Vec4 vecRecoil, vecRecoilRest;  // four-vectors of recoil momentum in fluid rest frame and lab frame
    Vec4 vecThermal, vecThermalLab;
    Parton newOne, newHole;         // parton object for the additionally created parton
    Source newSource;
    
    int imax = plist->size();       // number of partons in the list
    
    int id, newID;                  // original and new parton's ID
    int col, acol, newCol, newAcol; // number of the color string attached to the old, and new parton
    int ix, iy, iz, itau;           // parton's position in cell coordinates
    int ixmax, izmax;               // maximum cell number in x- and y-direction, and z-direction
    int radiate;                    // switches to hold what process will happen
    int radiatePhoton;
    int radiate_ggg;
    int radiate_gqq;
    int elastic_qq;
    int elastic_qg;
    int elastic_gq;
    int elastic_gg;
    int convert_quark;
    int convert_quark_to_gamma;
    int convert_gluon;
    int GluonToQuark;              // holds the decision whether a gluon converts to a quark (=1) or anti-quark (=0)
    //Sangyong's addition
    int mother;
    int daughter;
    int temp_pos;
    int in_coherence;              // flag for whether mother and daughter are still together
    int ind_proc;                  // which process to take place
    
    double dt = dtfm/hbarc;        // time step in [GeV^(-1)] - hbarc is defined in Basics.h 
                                   // Gamma = \int dGamma dp is given in units of GeV, so dt*Gamma is dimensionless
    double T;                      // temperature
    double p, q, qt;               // parton momentum in [GeV]
    double x, y, z;                // parton's position in [fm]
    double px, py, pz;             // parton's momentum
    double t, tau, eta;            // lab time and lab tau
    double vx, vy, vz, veta;       // flow velocity of the cell
    double vetaZ, vetaZTau;        // flow velocity of the cell
    double beta;                   // absolute value of the flow velocity
    double gamma;                  // gamma factor 
    double cosPhi, cosPhiRest;     // angle between pvec and v, and angle between pvec_restframe and v
    double cosPhiRestEl;           // angle between pvec_restframe and v in case of elastic collisions with transv. mom. transfer
    double boostBack;              // boost factor to boost back from rest- to lab-frame
    double pRest;                  // rest frame abs. value of momentum
    double pxRest, pyRest, pzRest; // rest frame momentum components
    double pRecoil;                // lab frame abs. value of recoil particle
    double newpx, newpy, newpz;    // temporary momentum components
    double pxRecoil, pyRecoil, pzRecoil;
    double pxThermalLab, pyThermalLab, pzThermalLab;
    double omega;                  // energy transfered in elastic collision
    double totalQuarkProb = 0.;    // total probability that the quark undergoes some interaction
    double totalGluonProb = 0.;    // total probability that the gluon undergoes some interaction
    double randx;                  // catch the random number
    double EMProb;                 // temporary holder for photon radiation rate
    double AccProb;                // accumulated probability
    double NextProb;               // next probability to consider
    double delp, k;                // sampled momentum in splitting and momentum of new parton
    double pdummy = 1.;
    double dE, dpx, dpy, dpz;      // energy deposited to medium
    
    ixmax = static_cast<int>(2.*hydroXmax/hydroDx+0.0001);
    izmax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);
    
    
    for ( int i=0; i<imax; i++)    // loop over all partons
    { 
        elastic_qq = 0;            // in the beginning assume nothing will happen
        elastic_qg = 0;            // then decide according to rates if processes occur 
        elastic_gq = 0;            // then decide according to rates if processes occur 
        elastic_gg = 0;
        radiate = 0;
        radiatePhoton = 0;
        radiate_ggg = 0;
        radiate_gqq = 0;
        convert_quark_to_gamma = 0;
        convert_quark = 0;
        convert_gluon = 0;
        
        counter++;   // counter that characterizes every step uniquely
                     // it is used to give new partons a unique color index
        if (plist->at(i).frozen() == 2) 
            continue; // if a particle is frozen out, don't evolve it any further
        
        id = plist->at(i).id();   // id of parton i 
        vecp = plist->at(i).p();  // four-vector p of parton i
        p = vecp.pAbs();          // |p| of parton i 
        
        x = plist->at(i).x();     // x value of position in [fm]
        y = plist->at(i).y();     // y value of position in [fm]
        z = plist->at(i).z();     // z value of position in [fm]
        
        // boost - using original position...
        ix = floor((hydroXmax+x)/hydroDx+0.0001);    // x-coordinate of the cell we are in now
        iy = floor((hydroXmax+y)/hydroDx+0.0001);    // y-coordinate of the cell we are in now
                                                     // note that x and y run from -hydroXmax to +hydroXmax
                                                     // and ix and iy from 0 to 2*hydroXmax hence the (hydroXmax+x or y) for both
        t = it*dtfm; // now let partons travel doing their vacuum shower stuff first...
        
        double pz_local = vecp.pz();
        double rapidity = 0.5*log( (p+pz_local)/(p-pz_local) );
        if(trackHistory and abs(rapidity)<1.0)
        {
            Tt->push_back(t);
            Tx->push_back(x);
            Ty->push_back(y);
            Tz->push_back(z);
            TE->push_back(p);
            Tpx->push_back(vecp.px());
            Tpy->push_back(vecp.py());
        }
        
        double eps = 1e-15;
        if (tauEtaCoordinates == 1) 
        {
            eta = 0.;
            if (fabs(t) < eps and fabs(z) < eps) 
                eta = 0.;
            else if ( t > 0 and fabs(z)<t ) 
                eta = 0.5 * log ( (t+z)/(t-z) );    //z is converted to eta in tau-eta coordinates
            else 
            {
                plist->at(i).frozen(2);
                continue;
            }
        }
        iz = floor((hydroZmax+eta)/hydroDz+0.0001); // z-coordinate of the cell we are in now
        // hydroZmax, hydroDz are actually hyDroEtamax, hydroDeta in tau-eta coordinates

        if (fixedTemperature == 0)
        {
            tau = 0.;           //just give tau some value - will change below
            plist->at(i).tFinal(t); // update the time - finally will be the final time
            plist->at(i).z(z+vecp.pz()/vecp.pAbs()*dtfm);      // update z position "
            // +++ if you want to move the partons before tau0 do the position update here.
            if ( moveBeforeTau0 == 1 )
            {
                plist->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
                plist->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
            }
            
            if(z*z<t*t)
                tau = sqrt(t*t-z*z);    // determine tau (z*z) should always be less than (t*t)
            else 
                continue;
            
            if (tau < hydroTau0)
            {
                continue;       // do not evolve if tau < tau0 or tau > tauMax
            }
            if (tau > hydroTauMax - hydroDtau)
            {
                plist->at(i).frozen(2);
                continue;       // do not evolve if tau < tau0 or tau > tauMax
            }
            
            // +++ if you DO NOT want to move the partons before tau0 do the position update here.
            if ( moveBeforeTau0 == 0 )
            {
                plist->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
                plist->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
            }
            
            // stop evolution if the position of the parton isout of range
            if ( ix < 0 || ix >= ixmax || iy < 0 || iy >= ixmax ||  iz < 0 || iz >= izmax ) 
            {
                plist->at(i).frozen(2);
                continue;
            }
            
            // get the temperature and flow velocity at the current position and time:
            hydroInfo = hydroSetup->getHydroValues(x, y, z, t);
            
            T = hydroInfo.T;
            vx = hydroInfo.vx;
            vy = hydroInfo.vy;
            vz = hydroInfo.vz;
            
            // parton stops interacting with medium if T < hydroTfinal
            // but still evolves
            if ( T < hydroTfinal ) continue;
            
            // the initial momentum vector
            px = vecp.px();
            py = vecp.py();
            pz = vecp.pz();
            
            // absolute value of flow velocity
            beta = sqrt(vx*vx+vy*vy+vz*vz);
            // angle between the particle and flow directions.
            cosPhi = (px*vx + py*vy + pz*vz)/(vecp.pAbs()*beta);
            gamma = 1./sqrt(1.-beta*beta);
            // boost p to fluid cell's rest frame
            pRest = p * gamma * (1.-beta*cosPhi);

            if(pCutPropToT)
                pCut = eLossCut*T; //standard: 4*T
            
            double eLossCut0 = eLossCut;
            // special treatment for negative particles
            if(plist->at(i).recoil() == -1) 
                eLossCut0 = 1.5;
            
            // parton stops evolution if energy is below threshold
            if(pRest/T < eLossCut0) 
                continue;
            
            pxRest = -vx*gamma*p 
              + (1.+(gamma-1.)*vx*vx/(beta*beta))*px 
              + (gamma-1.)*vx*vy/(beta*beta)*py
              + (gamma-1.)*vx*vz/(beta*beta)*pz;
            pyRest = -vy*gamma*p 
              + (1.+(gamma-1.)*vy*vy/(beta*beta))*py 
              + (gamma-1.)*vx*vy/(beta*beta)*px
              + (gamma-1.)*vy*vz/(beta*beta)*pz;
            pzRest = -vz*gamma*p 
              + (1.+(gamma-1.)*vz*vz/(beta*beta))*pz 
              + (gamma-1.)*vx*vz/(beta*beta)*px
              + (gamma-1.)*vy*vz/(beta*beta)*py;
            
            vecpRest.px(pxRest);
            vecpRest.py(pyRest);
            vecpRest.pz(pzRest);
            vecpRest.e(sqrt(pxRest*pxRest + pyRest*pyRest + pzRest*pzRest));
            
            // angle between pRest and flow vel.
            cosPhiRest = (pxRest*vx + pyRest*vy + pzRest*vz)/(pRest*beta);
        }
        else // the fixed temperature case:
        {
            plist->at(i).tFinal(t); // update the time - finally will be the final time
            beta = 0.;
            gamma = 1.; 
            cosPhiRest=1.;
            pRest = p;
            T = fixedT;
            vecpRest=vecp;
            plist->at(i).x(x+vecp.px()/vecp.pAbs()*dtfm);      // update x position for a massless parton (vel=c)
            plist->at(i).y(y+vecp.py()/vecp.pAbs()*dtfm);      // update y position "
            plist->at(i).z(z+vecp.pz()/vecp.pAbs()*dtfm);      // update z position "
            x = plist->at(i).x();                              // x value of position in [fm]
            y = plist->at(i).y();                              // y value of position in [fm]
            z = plist->at(i).z();                              // z value of position in [fm]
            if (pRest < eLossCut) // momentum cut: to avoid sampling issues later on  RMY Jan 11 2022
            {
                plist->at(i).frozen(2);
                continue;
            }
        }
      
        boostBack = gamma*(1.+beta*cosPhiRest);  // boost factor back from rest- to lab-frame 
        
        // Sangyong's addition starts
        // see if parton i has a mother or a daughter
        in_coherence = 0; // default no
        if (formationTime)
        {
            if(plist->at(i).daughter() != -1) 
                in_coherence = hasDaughter(plist, i, T);
            if(plist->at(i).mother() != -1) 
                in_coherence = hasMother(plist, i, T);
        }
        // has... returns 1 if still coherent
        // has... returns 0 if separated. also unhooks them.
        // Sangyong's addition ends

        alpha_s_rad = alpha_s;
        alpha_s_elas = alpha_s;
        
        if(runningRad || runningElas)
        {
            double alpha_s_max = 0.42;
            double alpha_s_min = 0.15;
            double beta0 = 9./(4.*M_PI);
            double lambda_QCD = 0.20;
            double eps = 1e-1;
            double Ca = 0.;
            if (fabs(id) > 0 and fabs(id) < 4) Ca = 4./3.;
            else if (id == 21) Ca = 3.;
            double mD_sq = (2./3.)*M_PI*alpha_s*pow(T, 2.)*(2.*Nc+Nf);
            double qmax_sq = 6.*pRest*T;
            double q_hat = Ca*alpha_s*mD_sq*T*log(1+qmax_sq/mD_sq);
        
            if(runningRad)
            {
                double mu = renormFacRad*pow(q_hat*pRest, 0.25);
                if(mu < lambda_QCD) mu = lambda_QCD + eps;
                alpha_s_rad = 1./(beta0*log(pow(mu/lambda_QCD, 2.)));
                
                if(alpha_s_rad < 0.)
                {
                  cout << "[Warning]::alpha_s_rad is smaller than 0." << endl;
                  alpha_s_rad = 0.;
                }
                if(alpha_s_rad > alpha_s_max) alpha_s_rad = alpha_s_max;
            }
        
            if(runningElas)
            {
                double qmin_sq =0.05*T;
                double scat_rate = Ca*alpha_s*T*(log(qmax_sq/qmin_sq) + 
                                                 log((qmin_sq+mD_sq)/(qmax_sq+mD_sq)));
                double mfp = 1./scat_rate;
                double mu = renormFacElas*sqrt(q_hat*mfp);
                if(mu < lambda_QCD) mu = lambda_QCD + eps;
                alpha_s_elas = 1./(beta0*log(pow(mu/lambda_QCD, 2.)));
              
                if(alpha_s_elas < 0.)
                {
                  cout << "[Warning]::alpha_s_elas is smaller than 0." << endl;
                  alpha_s_elas = 0.;
                }
                if(alpha_s_elas > alpha_s_max)       alpha_s_elas = alpha_s_max;
                else if (alpha_s_elas < alpha_s_min) alpha_s_elas = alpha_s_min;
            }
        }

        // get total probabilities for radiative processes
        if (doRadiative == 1 and pRest/T > AMYpCut) 
        {
            n = rates->integrate(pRest/T,T,alpha_s_rad, Nf);
            // n stores the areas under the prob. distr.
            // Warn if probability for any process is larger than 1:
            if (dt*n.Gamma/gamma>1. and pRest/T>AMYpCut) 
                cout << fixed << "WARNING: probability to emit during one time step " << dt*n.Gamma/gamma 
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*n.Gamma_ggg/gamma >1. and pRest/T>AMYpCut) 
                cout << fixed << "WARNING: probability ggg during one time step " << dt*n.Gamma_ggg/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << ", T=" << T << endl;
            if (dt*n.Gamma_gqq/gamma>1. and pRest/T>AMYpCut) 
                cout << fixed << "WARNING: probability gqq during one time step " << dt*n.Gamma_gqq/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*n.Gamma_em/gamma>1. and pRest/T>AMYpCut) 
                cout << fixed << "WARNING: probability qqgamma during one time step " << dt*n.Gamma_em/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
        }
      
        // get total probabilities for elastic processes

        if (doElastic == 1 ) 
        {
            qqRate = elastic->totalRate(pRest, T, alpha_s_elas, Nf, 1);
            gqRate = qqRate*9./4.;
            qgRate = elastic->totalRate(pRest, T, alpha_s_elas, Nf, 3);
            ggRate = qgRate*9./4.; 
            conversionqg = rates->Gammaqg( pRest, T, alpha_s_elas );
            conversiongq = rates->Gammagq( pRest, T, alpha_s_elas, Nf );
            if ( photonSwitch == 1 ) 
            {
                conversionqgamma = rates->Gammaqgamma( pRest, T, alpha_s );
                if ( abs(id) == 2 ) 
                    conversionqgamma*=(4./9.); // multiplying by (ef/e)^2
                else if ( abs(id) == 1 or abs(id) == 3 ) 
                    conversionqgamma*=(1./9.); // multiplying by (ef/e)^2
            }
            else 
                conversionqgamma = 0.;
        
            // Warn if probability for any process is larger than 1:
            if (dt*qqRate/gamma>1.) 
                cout << "WARNING: probability qq during one time step " << dt*qqRate/gamma 
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*gqRate/gamma >1.) 
                cout << "WARNING: probability gq during one time step " << dt*gqRate/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*qgRate/gamma>1.) 
                cout << "WARNING: probability qg during one time step " << dt*qgRate/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*ggRate/gamma>1.) 
                cout << "WARNING: probability gg during one time step " << dt*ggRate/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*conversionqg/gamma>1.) 
                cout << "WARNING: probability conversion q->g during one time step " << dt*conversionqg/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*conversionqgamma/gamma>1.) 
                cout << "WARNING: probability conversion q->gamma during one time step " << dt*conversionqgamma/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
            if (dt*conversiongq/gamma>1.) 
                cout << "WARNING: probability conversion g->q during one time step " << dt*conversiongq/gamma
                     << " > 1. Decrease time step to improve result. pRest/T=" << pRest/T << endl;
        }
  

      if ( abs(id) > 0 and abs(id) < 4 ) // if parton is a quark (u, d, or s), let it evolve like a quark
      {
        // Check if we want a perturbative calculation. If so, do that here
        // RY June 11 2021
        if ( pertCalc ) 
        {
            // Commented out: Just create the photons and add them to the parton list
            // negative: the parton list can become too large, positive: more control at the 
            // analysis stage. Rouz, August 2021 
            // first ensure that both histograms for perturbative calculation
            // have been allocated. They should be.
            //if ( conversion_photons != NULL and amy_photons != NULL)
            //{
                
                //Parton part = plist->at(i);
                //double chargeSquared = abs(id) == 2 ? 4/9. : 1/9.;
                //double boostFactor =  gamma*(1.0 - beta*cosPhi);
                //double weightMartini;//, weightAngular;
          
                //if (part.status() > 0 ) //&& useRateTableforConversion)
                //{
                //    if(abs(part.p().eta()) < etaCut_pert)
                //    // if the particle is in the pseudorapidity window, apply collinear approximation
                //    {
                //        weightMartini = chargeSquared*rates->Gammaqgamma( pRest, T, alpha_s );
                //        conversion_photons->Fill(part.p().pT(), weightMartini*boostFactor*dt/gamma);
                //    }
                //    if(pRest/T > AMYpCut)
                //    {   
                //        Norms tempNorm = rates->integrate(pRest/T, T, alpha_s, Nf);
                //        ReturnValue dd = rates->findValuesRejection(pRest/T, T, alpha_s, Nf, random, import, 4);
                //        Parton pp;
                //        pp.p(dd.y*T*vecp.px()/vecp.pAbs()*boostBack,  // emitted gamma's momentum
                //             dd.y*T*vecp.py()/vecp.pAbs()*boostBack,  // since the direction does not change I can use vecp here already
                //             dd.y*T*vecp.pz()/vecp.pAbs()*boostBack);
                //        if ( abs(pp.p().eta()) < etaCut_pert)
                //        {
                //             amy_photons->Fill(pp.p().pT(), chargeSquared*boostFactor*tempNorm.Gamma_em*dt/gamma);
                //        } 
                //    }
                //}
            //}
            //else
            //{
            //    cout<<"ERROR: Histograms for perturbative calculation are not allocated. EXIT."<<endl;
            //    exit(-1);
            //}
            //Deal with the bremsstrahlung photon first:
            int n_photons;
            if (pRest/T > AMYpCut)
            {            
                double brem_rate = rates->integrate( pRest, T, alpha_s_rad, Nf).Gamma_em*exaggerate_brem*dt/gamma;
                brem_rate = abs(id) == 2 ? brem_rate*(4./9) : brem_rate*(1./9);
                double brem_rate_compare = random->genrand64_real1();
                n_photons = floor(brem_rate);
                if (brem_rate - n_photons > brem_rate_compare)
                    n_photons += 1;
                ReturnValue d;
                double DELTAP;
                for (int j = 0; j < n_photons; j++)
                { 
                    d = rates->findValuesRejection(pRest/T, T, alpha_s_rad, Nf, random, import, 4); 
                    DELTAP = d.y*T;
                    //cout << "\t brem photon number "<< j << " delp: "<< delp << endl;
                    if ( DELTAP < pRest and DELTAP > 0. and DELTAP > pCut)
                    {
                        //p = (pRest - DELTAP)*boostBack;
                        //k = DELTAP;
                        Parton photonAMY;
                        photonAMY.p(DELTAP*vecp.px()/vecp.pAbs()*boostBack,  // emitted gamma's momentum
                                    DELTAP*vecp.py()/vecp.pAbs()*boostBack,
                                    DELTAP*vecp.pz()/vecp.pAbs()*boostBack);
                        photonAMY.id(22);        // emitted parton is a photon
                        photonAMY.mass(0.);
                        photonAMY.frozen(2);     // photons do not interact anymore
                        photonAMY.x(x);          // set the new parton's initial position
                        photonAMY.y(y);
                        photonAMY.z(z);
                        photonAMY.tFinal(t);    // creation time of this photon
                        photonAMY.tauAtEmission(tau); //creation (proper) time of the photon
                        photonAMY.source(2);     // AMY photon
                        photonAMY.daughter_of(-1); //don't update mother and daughter information
                        photonAMY.mother_of(-1);  // this is a perturbative calculation, piggy-backing off the main MC engine
                        photonAMY.recoil(plist->at(i).recoil());//inherit the recoil code of mother
                        plist->push_back(photonAMY); // add the photon to the list of partons
                    }
                }
            }
            //cout<<"num brem photons : "<<n_photons<<endl;
            //Deal with conversion photon next:
            double convrate = rates->Gammaqgamma(pRest, T, alpha_s_elas)*exaggerate_conv*dt/gamma;
            convrate = abs(id) == 2 ? convrate*(4./9) : convrate*(1./9);
            double conv_rate_compare = random->genrand64_real1();
            n_photons = floor(convrate);
            if (convrate - n_photons > conv_rate_compare)
                n_photons += 1;
            //cout << "num conv photons : "<<n_photons<<endl;
            for ( int j = 0 ; j < n_photons; j++)
            {
                Parton photonConv;
                photonConv.col(0);
                photonConv.acol(0);
                photonConv.id(22);
                photonConv.mass(0.);
                photonConv.frozen(2); 
                photonConv.x(x);
                photonConv.y(y);
                photonConv.z(z);
                photonConv.tFinal(t); // creation time of this photon
                photonConv.tauAtEmission(tau);
                photonConv.p(px, py, pz);
                photonConv.daughter_of(-1);
                photonConv.mother_of(-1);
                photonConv.source(1);//conversion photon
                photonConv.recoil(plist->at(i).recoil()); // inherit the recoil code of mother
                plist->push_back(photonConv);
            }
            // Done with the perturbative stuff.
        }    
          // in the following decide what to do in this time step
          // Here begin the modifications that will allow elastic collisions below the 
          // AMYpCut, while not sampling radiative splittings below this scale. -CFY 9/28/2011
          if (doRadiative == 1 and doElastic == 1 ) //radiative + elastic 
          {
              if ( abs(id) == 2 ) // multiply by (e_f/e)^2
                  n.Gamma_em*=4./9.;
              else
                  n.Gamma_em*=1./9.;

              if(pRest/T>AMYpCut)
              {
                  totalQuarkProb = dt/gamma*(n.Gamma)
                                 + dt/gamma*(qqRate)+dt/gamma*(qgRate)
                                 + dt/gamma*conversionqg+dt/gamma*conversionqgamma;
              }
              else
              {
                  totalQuarkProb = dt/gamma*(qqRate)+dt/gamma*(qgRate)
                                 + dt/gamma*conversionqg+dt/gamma*conversionqgamma;
              }

              if ( photonSwitch == 1 && pRest/T>AMYpCut) 
              {
                  totalQuarkProb += dt/gamma*(n.Gamma_em);
                  EMProb = (dt/gamma*(n.Gamma_em))/totalQuarkProb;
              }
              else
              {
                  EMProb = 0.0;
              }

              if (totalQuarkProb>1)
              {
                  cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
                       << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
                  cout << "Parton with pRest : "<<pRest<<" and source : "<<plist->at(i).source()<<endl;
              }
        
              if ( random->genrand64_real1() < totalQuarkProb ) // check if something happens with the quark
              {
                  
                  if(pRest/T>AMYpCut)
                  {
                      // Sangyong's addition:
                      // in_coherence is the new flag
                      // added if(in_coherence==0) to radiate* and convert*
                      // this includes coherence effects in radiations, but not in conversions
                      // aways do elastic - elastic was lumped together with the conversions
                              // I am modifying this to use one random number
                      // we either have
                              // totalQuarkProb = dt/gamma*(n.Gamma)
                      //  +dt/gamma*(qqRate)+dt/gamma*(qgRate)
                      //  +dt/gamma*conversionqg+dt/gamma*conversionqgamma; 
                      // or
                              // totalQuarkProb = dt/gamma*(n.Gamma)
                      //  +dt/gamma*(qqRate)+dt/gamma*(qgRate)
                      //  +dt/gamma*conversionqg+dt/gamma*conversionqgamma; 
                      //  +dt/gamma*(n.Gamma_em);
                      
                      randx = random->genrand64_real1();
                      
                      AccProb = 0.0;
                      NextProb = (dt/gamma*(qqRate))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                          elastic_qq = 1;

                      AccProb += NextProb;
                      NextProb = (dt/gamma*(qgRate))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                          elastic_qg = 1;
                      
                      AccProb += NextProb;
                      NextProb = dt/gamma*(n.Gamma)/totalQuarkProb;
                      if(AccProb <= randx && randx < (AccProb + NextProb))
                      {
                          if(in_coherence==0) 
                              radiate = 1;
                      }

                      AccProb += NextProb;
                      NextProb = (dt/gamma*conversionqg)/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                      {
                          if(in_coherence==0 or plist->at(i).daughter() != -1)
                              convert_quark = 1;
                      }

                      // +dt/gamma*conversionqgamma; 
                      AccProb += NextProb;
                      NextProb = (dt/gamma*conversionqgamma)/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                      {
                          if(in_coherence==0 or plist->at(i).daughter() != -1)
                             convert_quark_to_gamma = 1;
                      }

                      // EMProb = dt/gamma*(n.Gamma_em)/totalQuarkProb
                      // this is set to zero when photonSwitch == 0
                      AccProb += NextProb;
                      NextProb = EMProb;
                      if( AccProb <= randx && randx <= (AccProb + NextProb) )
                      {
                          if(in_coherence==0) 
                              radiatePhoton = 1;
                      }
                  }// if(pRest/T>AMYpCut)
                  else
                  {
                  randx = random->genrand64_real1();

                      AccProb = 0.0;
                      NextProb = (dt/gamma*(qqRate))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                          elastic_qq = 1;

                      AccProb += NextProb;
                      NextProb = (dt/gamma*(qgRate))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                          elastic_qg = 1;

                      AccProb += NextProb;
                      NextProb = (dt/gamma*(conversionqg))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                      {
                          if(in_coherence==0 or plist->at(i).daughter() != -1)
                              convert_quark = 1;
                      }

                      AccProb += NextProb;
                      NextProb = (dt/gamma*(conversionqgamma))/totalQuarkProb; 
                      if( AccProb <= randx && randx < (AccProb + NextProb) )
                      {
                          if(in_coherence==0 or plist->at(i).daughter() != -1)
                              convert_quark_to_gamma = 1;
                      }
                  }// pRest/T < AMYpCut
              }// something happens to the quark
              else  // nothing happens to the quark
              {
                  continue;
              }// nothing happens
          }// rad + elastic
          //Here ends the modification of the rates for quarks, on to gluons. -CFY
          else if (doRadiative == 1 and doElastic == 0 and pRest/T>AMYpCut) //radiative only
          {
              // Sangyong's addition:
              // in_coherence is the new flag defined in Parton.h
              // originally, I added if(in_coherence==0) to radiate*
              // but since this part is radiative only, that does not make sense. removed.
              totalQuarkProb = dt/gamma*(n.Gamma);
              if ( photonSwitch == 1 ) totalQuarkProb += dt/gamma*(n.Gamma_em);

              if (totalQuarkProb>1)
                  cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
                       << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
              if ( random->genrand64_real1() < totalQuarkProb ) // check if something happens with the quark
              {
                  if ( random->genrand64_real1() < dt/gamma*(n.Gamma)/totalQuarkProb )
                      radiate = 1; 
                  else if ( photonSwitch == 1 )
                      radiatePhoton = 1;
              }
              else // nothing happens to the quark
              {
                  continue;
              }
          } // rad only
          else if (doRadiative == 0 and doElastic == 1) //elastic only
          // no sense in implementing coherence effect here, either
          {
                totalQuarkProb = dt/gamma*(qqRate)+dt/gamma*(qgRate)
                   + dt/gamma*conversionqg+dt/gamma*conversionqgamma; 
                // dt/gamma is dt_rest-frame
                
                if (totalQuarkProb>1)
                cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
                     << totalQuarkProb << " > 1. Reduce time step to cure this." << endl; 
                if ( random->genrand64_real1() < totalQuarkProb )  // check if something happens
                {
                    if ( random->genrand64_real1() < conversionqg/(qqRate+qgRate+conversionqg+conversionqgamma) )
                    {
                        convert_quark = 1;
                    }
                    else if ( random->genrand64_real1() < conversionqgamma/(qqRate+qgRate+conversionqgamma) )
                    {
                        convert_quark_to_gamma = 1;
                    }
                    else 
                    {
                        if ( random->genrand64_real1() < qqRate/(qqRate+qgRate) )
                        {
                            elastic_qq = 1;
                            elastic_qg = 0;
                        }
                        else
                        {
                            elastic_qq = 0;
                            elastic_qg = 1;
                        }
                    }
                }
                else 
                {
                    continue;
                }
          }
          else
              continue;

          // now do what has been decided before
          if( radiate == 1 ) // see if emission happens
          {
              if(pRest/T>AMYpCut) // do not evolve partons with momenta below this scale 
              {
                  // do process 1, q->qg
                  f = rates->findValuesRejection(pRest/T, T, alpha_s_rad, Nf, random, import, 1); 
                  delp = f.y*T;
                  if(delp > pRest) continue;

                  p = (pRest - delp)*boostBack;    // quark's new momentum in the lab frame
                  newvecp=p/vecp.pAbs()*vecp;      // new quark momentum
                  plist->at(i).p(newvecp);         // change quark's momentum to the new one

                  if(delp > pRest - pCut)
                  {
                      plist->at(i).frozen(2);
                      // exclude soft recoil by assigning id == 0
                      if(plist->at(i).recoil() == -1) plist->at(i).id(0);
                  }

                  if (delp > pCut) // if new parton is kept (f.y*T > threshold [in GeV])
                  {
                      newOne.p(delp*vecp.px()/vecp.pAbs()*boostBack,  // emitted gluon's momentum
                               delp*vecp.py()/vecp.pAbs()*boostBack,  // since the direction does not change I can use vecp here already
                               delp*vecp.pz()/vecp.pAbs()*boostBack);
                      col = plist->at(i).col();              // color of original quark  
                      acol = plist->at(i).acol();            // anti-color
                      plist->at(i).increment_hard_radiation(); // Rouz : count hard radiations
                      if (col!=0)                               // if we had a quark
                      {
                          plist->at(i).col(100000000+counter);  // set color to new color
                          newOne.col(col);                      // set new gluon's color to quark's original color
                          newOne.acol(100000000+counter);       // set new gluon's anti-color to quark's new color
                      }
                      else if (acol!=0) // if we had an anti-quark
                      {
                          plist->at(i).acol(100000000+counter);  // set anti-color to new color
                          newOne.col(100000000+counter); // set new particle's color
                          newOne.acol(acol);            // set new particle's anti-color   
                      }

                      newOne.id(21);    // emitted parton is a gluon
                      newOne.init_counts();//Rouz : initialize the counts
                      newOne.mass(0.);
                      newOne.frozen(0);
                      newOne.x(x);      // set the new parton's initial position
                      newOne.y(y);
                      newOne.z(z);
                      newOne.source(11); // AMY parton, gluon emitted from quark or anti-quark
                      newOne.recoil(plist->at(i).recoil());
          
                      //Sangyong's addition for splits
                      daughter = plist->size(); // the position of daughter in the plist when push_back'ed
                      newOne.daughter_of(i); // daughter of the current parton
                      newOne.mother_of(-1);  // mother of nobody
                      newOne.p_at_split(newOne.p()); // momentum at the split point
                      
                      plist->at(i).daughter_of(-1); // daughter of nobody
                      plist->at(i).mother_of(daughter);
                      plist->at(i).p_at_split(plist->at(i).p());
                      //Sangyong's addition end
                      plist->push_back(newOne);  // add the gluon to the list of partons if (f.y>AMYpCut)
                  }
                  else
                  {
                    // delp < pCut: soft radiation.
                    plist->at(i).increment_soft_radiation();
                  }
                  if(outputSource)
                  {
                      // if delp < 0 or delp > pRest, jet energy gain from thermal medium.
                      // dE is negative but the direction is same as mother jet.
                      // if 0 < delp < pCut or pRest - pCut < delp < pRest, 
                      // radiated jet particle is absorbed into thermal medium.
                      if(delp < pCut)
                      {
                          dE = delp;
                          dpx = fabs(delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(delp)/vecp.pAbs()*vecp.pz();
                      }
                      else if (delp > pRest - pCut)
                      {
                          dE = pRest-delp;
                          dpx = fabs(pRest-delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(pRest-delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(pRest-delp)/vecp.pAbs()*vecp.pz();
                      }
                      newSource.type(2);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
              }
          }
          else if( radiatePhoton == 1 ) // see if photon emission happens
          {
              if(pRest/T>0.01) // do not evolve partons with momenta below this scale
              {
                  // do process 4, q->qgamma
                  f = rates->findValuesRejection(pRest/T, T, alpha_s_rad, Nf, random, import, 4); 
                  delp = f.y*T;

                  if(delp > pRest || delp < 0.) continue;           // do nothing at this moment
                  else if(delp > pRest - pCut)
                  {
                      plist->at(i).frozen(2);
                      // exclude soft recoil by assigning id == 0
                      if(plist->at(i).recoil() == -1) plist->at(i).id(0);
                  }
                  p = (pRest - delp)*boostBack;
                  k = delp;
                  newvecp = p/vecp.pAbs()*vecp;               // four vector for new quark momentum
                  plist->at(i).p(newvecp);                 // update new quark momentum
                   
                  newOne.p(k*vecp.px()/vecp.pAbs()*boostBack,  // emitted gamma's momentum
                           k*vecp.py()/vecp.pAbs()*boostBack,
                           k*vecp.pz()/vecp.pAbs()*boostBack);
                  newOne.id(22);                 // emitted parton is a photon
                  newOne.init_counts(); // Rouz's addition. New parton has had no scatterings so far
                  //newOne.splits(0);
                  newOne.mass(0.);
                  newOne.frozen(2);              // photons do not interact anymore
                  newOne.x(x);                   // set the new parton's initial position
                  newOne.y(y);
                  newOne.z(z);
                  newOne.tauAtEmission(tau);
                  newOne.source(2);              // AMY photon
                  newOne.recoil(plist->at(i).recoil());

                  //Sangyong's addition for splits
                  daughter = plist->size(); // the position of daughter in the plist when push_back'ed
                  newOne.daughter_of(i); // daughter of the current parton
                  newOne.mother_of(-1);  // mother of nobody
                  newOne.p_at_split(newOne.p()); // momentum at the split point
                  
                  plist->at(i).daughter_of(-1); // daughter of nobody
                  plist->at(i).mother_of(daughter);
                  plist->at(i).p_at_split(plist->at(i).p());
                  //Sangyong's addition

                  plist->push_back(newOne); // add the photon to the list of partons
                  if (delp < pCut) // Rouz : update count of soft or hard radiations
                    plist->at(i).increment_soft_radiation();
                  else
                    plist->at(i).increment_hard_radiation();

                  if(outputSource and delp > pRest - pCut)
                  {
                      // if delp > pRest, jet is converted to photon
                      // and gains energy from thermal medium
                      // dE is negative but the direction is same as mother jet.
                      // if pRest - pCut < delp < pRest, 
                      // radiated jet particle is absorbed into thermal medium.
                      dE = pRest-delp;
                      dpx = fabs(pRest-delp)/vecp.pAbs()*vecp.px();
                      dpy = fabs(pRest-delp)/vecp.pAbs()*vecp.py();
                      dpz = fabs(pRest-delp)/vecp.pAbs()*vecp.pz();

                      newSource.type(2);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);
                      slist->push_back(newSource);
                  }
              }
          }                             
          else if ( convert_quark == 1  ) // do q->g
          {
              col = plist->at(i).col();   // color of original quark  
              acol = plist->at(i).acol(); // anti-color
              if (col!=0)                 // if we had a quark
              {
                  plist->at(i).acol(200000000+counter);  // add a new anti-color for the gluon
                  newOne.col(200000000+counter);         // create new thermal gluon
                  newOne.acol(300000000+counter);                        
              }
              else if (acol!=0)
              {
                  plist->at(i).col(200000000+counter);   // add a new color for the gluon
                  newOne.acol(200000000+counter);        // create new thermal gluon
                  newOne.col(300000000+counter);                        
              }
              plist->at(i).id(21);    // convert to gluon
              plist->at(i).source(14); // Conversion gluon
              double choice = random->genrand64_real1();
              if ( choice < 0.5 )
              {
                newOne.id(21);          // add a thermal gluon
                newOne.mass(0.); 
                newOne.frozen(2);
                newOne.x(x);            // set the new parton's initial position
                newOne.y(y);
                newOne.z(z);
                newOne.p(random->thermal(T, -1));  // sample momentum from Bose distribution at temperature T
                newOne.recoil(plist->at(i).recoil());

                //Sangyong's addition for formation time of radiation
                newOne.daughter_of(-1); // daughter of nobody
                newOne.mother_of(-1);  // mother of nobody
                newOne.p_at_split(newOne.p()); // momentum at the conversion point

                plist->push_back(newOne); 
              }
              else 
              {
                 newOne.id(id);                        // add a thermal quark of same flavor
                 if (col!=0)                               // if we had a quark
                 {
                     newOne.col(300000000+counter);            // attach the color string to the thermal quark
                     newOne.acol(0);                        
                 }
                 else if (acol!=0)                         // if we had an anti-quark
                 {
                     newOne.acol(300000000+counter);           // attach the color string to the thermal anti-quark
                     newOne.col(0);                        
                 }
                 if(abs(id) < 3) 
                     newOne.mass(0.33);
                 else 
                     newOne.mass(0.5);

                 newOne.frozen(2);
                 newOne.x(x);                              // set the new parton's initial position
                 newOne.y(y);
                 newOne.z(z);
                 newOne.p(random->thermal(T, 1));          // sample momentum from Fermi distribution at temperature T
                 newOne.recoil(plist->at(i).recoil());
                 
                 //Sangyong's addition for formation time of radiation
                 newOne.daughter_of(-1); // daughter of nobody
                 newOne.mother_of(-1);  // mother of nobody
                 newOne.p_at_split(newOne.p()); // momentum at the conversion point
                 
                 plist->push_back(newOne); 
              }
              continue;                                 //prevent the new gluon from interacting again in this time step 
          }
          else if ( convert_quark_to_gamma == 1 and pRest/T>0.01 ) // do q->gamma
          {
              col = plist->at(i).col();  // color of original quark  
              acol = plist->at(i).acol();// anti-color
              plist->at(i).acol(0);      // make color neutral photon
              plist->at(i).col(0);       // make color neutral photon
              if (col!=0)                // if we had a quark
              {
                  newOne.col(col);                      // attach the color string to the thermal gluon
                  newOne.acol(400000000+counter);                        
              }
              else if (acol!=0)
              {
                  newOne.acol(acol);                    // attach the color string to the thermal gluon
                  newOne.col(400000000+counter);                        
              }

              plist->at(i).id(22);    // convert to photon
              plist->at(i).source(1); // Conversion photon "source" code.
              plist->at(i).init_counts(); // Rouz: this sets all these counts to zero for the photon. they don't interact in MARTINI
              //TODO: here, why not a coin flip for thermal q/qbar vs gluon? Rouz Dec 8, 2021
              newOne.id(21);          // add a thermal gluon
              newOne.init_counts(); // Rouz: new parton, hasn't had any interactions yet
              newOne.mass(0.); 
              newOne.frozen(0);
              newOne.x(x);            // set the new parton's initial position
              newOne.y(y);
              newOne.z(z);
              newOne.tauAtEmission(tau);
              newOne.p(random->thermal(T, -1)); // sample momentum from Bose distribution at temperature T
              newOne.recoil(plist->at(i).recoil());

              //Sangyong's addition for formation time of radiation
              newOne.daughter_of(-1); // daughter of nobody
              newOne.mother_of(-1);  // mother of nobody
              newOne.p_at_split(newOne.p()); // momentum at the conversion point

              plist->push_back(newOne); 

              newOne.id(id);  // add a thermal quark of same flavor as initial quark
              if (col!=0)     // if we had a quark
              {
                  newOne.col(400000000+counter); // attach the color string to the thermal quark
                  newOne.acol(0);                        
              }
              else if (acol!=0) // if we had an anti-quark
              {
                  newOne.acol(400000000+counter); // attach the color string to the thermal anti-quark
                  newOne.col(0);                        
              }
              if ( abs(id) < 3) 
                  newOne.mass(0.33);
              else 
                  newOne.mass(0.5);

              newOne.frozen(0);
              newOne.x(x); // set the new parton's initial position
              newOne.y(y);
              newOne.z(z);
              newOne.p(random->thermal(T, 1)); // sample momentum from Fermi distribution at temperature T
              newOne.recoil(plist->at(i).recoil());

              //Sangyong's addition for formation time of radiation
              newOne.daughter_of(-1); // daughter of nobody
              newOne.mother_of(-1);  // mother of nobody
              newOne.p_at_split(newOne.p()); // momentum at the conversion point
              
              plist->push_back(newOne); 
              continue;                                 //prevent the new gluon from interacting again in this time step 
          }
          else if ( elastic_qq == 1 ) // do qq->qq
          {
              for(;;)
              {
                  omega = elastic->findValuesRejection(pRest/T, T, alpha_s_elas, Nf, random, import, 1); // this is omega/T
                  if(fabs(omega) < pRest/T) break;
              }
              
              if( transferTransverseMomentum == 1 and omega < 800. )//&& fabs(eta) < 4.
              {
                  int repeat = 1;
                  // k_min is a mimumum momenum of a thermal parton to be sample 
                  // for a recoil parton to be on-shell
                  double k_min;
                  do
                  {
                      q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s_elas, Nf, random, import, 1);
                      k_min = (q-omega)/2.;
                      if(k_min < 20.) repeat = 0;
                  } while (repeat);
              }
              else 
              {
                  q = omega;
              }
          }
          else if ( elastic_qg == 1 ) // do qg->qg
          {
              for(;;)
              {
                  omega = elastic->findValuesRejection(pRest/T, T, alpha_s_elas, Nf, random, import, 3); // this is omega/T
                  if(fabs(omega) < pRest/T) break;
              }

              if( transferTransverseMomentum == 1 && omega < 800. ) // && fabs(eta) < 4.)
              {
                  int repeat = 1;
                  // k_min is a mimumum momenum of a thermal parton to be sample 
                  // for a recoil parton to be on-shell
                  double k_min;
                  do
                  {
                      q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s_elas, Nf, random, import, 3);
                      k_min = (q-omega)/2.;
                      if(k_min < 20.) repeat = 0;
                  } while (repeat);
              }
              else 
              {
                  q = omega;
              }
          }
          if ( elastic_qq == 1 || elastic_qg == 1 ) // in both cases ( qq->qq and qg->qg ) change the momentum
          {
              newvecpRest=elastic->getNewMomentum(vecpRest, omega*T, q*T, random); // takes everything in GeV
              if (beta > 1e-15) // before: beta == 0, should avoid equality check for doubles
              {
                  newpx = vx*gamma*newvecpRest.pAbs() 
                    + (1.+(gamma-1.)*vx*vx/(beta*beta))*newvecpRest.px() 
                    + (gamma-1.)*vx*vy/(beta*beta)*newvecpRest.py()
                    + (gamma-1.)*vx*vz/(beta*beta)*newvecpRest.pz();
                  newpy = vy*gamma*newvecpRest.pAbs() 
                    + (1.+(gamma-1.)*vy*vy/(beta*beta))*newvecpRest.py()
                    + (gamma-1.)*vx*vy/(beta*beta)*newvecpRest.px()
                    + (gamma-1.)*vy*vz/(beta*beta)*newvecpRest.pz();
                  newpz = vz*gamma*newvecpRest.pAbs() 
                    + (1.+(gamma-1.)*vz*vz/(beta*beta))*newvecpRest.pz()
                    + (gamma-1.)*vx*vz/(beta*beta)*newvecpRest.px()
                    + (gamma-1.)*vy*vz/(beta*beta)*newvecpRest.py();
                  
                  // momentum four-vector in lab frame
                  newvecp.px(newpx);
                  newvecp.py(newpy);
                  newvecp.pz(newpz);
                  newvecp.e(sqrt(newpx*newpx+newpy*newpy+newpz*newpz));
              }
              else
                newvecp = newvecpRest;
              double newpRest = newvecpRest.pAbs();
              if(newpRest < pCut)
              {
                  if(outputSource)
                  {
                      // output energy deposition from jet particle absorption
                      dE = newvecpRest.pAbs();
                      dpx = newvecpRest.px();
                      dpy = newvecpRest.py();
                      dpz = newvecpRest.pz();

                      newSource.type(1);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
                  plist->at(i).frozen(2);
                  //newvecp *= pdummy/newvecp.pAbs();
                  // exclude soft recoil by assigning id == 0
                  if(plist->at(i).recoil() == -1) plist->at(i).id(0);
              }
    
              plist->at(i).p(newvecp);   // change quark's momentum to the new one

              vecq = vecpRest - newvecpRest; // momentum transfer q
              if (vecq.pAbs() > pCut) // Rouz : Count elastic scatterings
                plist->at(i).increment_hard_elastic();
              else
                plist->at(i).increment_soft_elastic();

              if(getRecoil and omega > 0. and plist->at(i).recoil() != -1)
              {
                  if( elastic_qq == 1 )   // quark recoil
                  {
                      // thermal quark to be recoiled
                      vecThermal = elastic->getThermalMomentum(vecq, T, 1, random);
                      vecRecoilRest = vecq + vecThermal;

                      pxRecoil = vx*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecRecoilRest.px() 
                        + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.py()
                        + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.pz();
                      pyRecoil = vy*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecRecoilRest.py()
                        + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.px()
                        + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.pz();
                      pzRecoil = vz*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecRecoilRest.pz()
                        + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.px()
                        + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.py();
                      
                      // momentum four-vector in lab frame
                      vecRecoil.px(pxRecoil);
                      vecRecoil.py(pyRecoil);
                      vecRecoil.pz(pzRecoil);
                      vecRecoil.e(sqrt(pxRecoil*pxRecoil +
                                       pyRecoil*pyRecoil +
                                       pzRecoil*pzRecoil));

                      if( vecRecoilRest.pAbs() > recoilCut )
                      {
                          double r = random->genrand64_real1();
                          // note that when using general N_f this has to be changed
                          if (r<0.33) newID=1;
                          else if (r<0.66) newID=2;
                          else newID=3;
                          double mass;
                          if (newID<3) mass=0.33;  // set the quark's mass
                          else mass = 0.5;

                          newOne.p(vecRecoil);
                          newOne.mass(mass);

                          double matter = random->genrand64_real1();
                          if(matter < 0.5)   // thermal parton is quark
                          {
                              newOne.id(newID);
                              newOne.col(500000000+counter);
                              newOne.acol(0);
                              newHole.id(newID);
                        newHole.col(600000000+counter);
                        newHole.acol(0);
                          }
                          else  // thermal parton is anti-quark
                          {
                              newOne.id(-newID);
                              newOne.col(0);
                              newOne.acol(500000000+counter);
                              newHole.id(-newID);
                        newHole.col(0);
                        newHole.acol(600000000+counter);
                          }

                          newOne.x(x);
                          newOne.y(y);
                          newOne.z(z);
                          newOne.frozen(0);
                          newOne.recoil(1);
                          newOne.init_counts(); // Rouz: initialize the counts
                          //Sangyong's addition for formation time of radiation
                          newOne.daughter_of(-1); // daughter of nobody
                          newOne.mother_of(-1);  // mother of nobody
                          newOne.p_at_split(newOne.p()); // momentum at the conversion point

                          plist->push_back(newOne);    // add recoil parton to the parton list

                          pxThermalLab = vx*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecThermal.px()
                                    + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.py()
                                    + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.pz();
                          pyThermalLab = vy*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecThermal.py()
                                    + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.px()
                                    + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.pz();
                          pzThermalLab = vz*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecThermal.pz()
                                    + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.px()
                                    + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.py();

                          // momentum four-vector in lab frame
                          vecThermalLab.px(pxThermalLab);
                          vecThermalLab.py(pyThermalLab);
                          vecThermalLab.pz(pzThermalLab);
                          vecThermalLab.e(sqrt(pxThermalLab*pxThermalLab +
                                                pyThermalLab*pyThermalLab +
                                                pzThermalLab*pzThermalLab));

                          newHole.p(vecThermalLab);
                          newHole.mass(mass);

                          newHole.x(x);
                          newHole.y(y);
                          newHole.z(z);
                          newHole.frozen(0);
                          newHole.recoil(-1);

                          //Sangyong's addition for formation time of radiation
                          newHole.daughter_of(-1); // daughter of nobody
                          newHole.mother_of(-1);  // mother of nobody
                          newHole.p_at_split(newHole.p()); // momentum at the conversion point

                          plist->push_back(newHole);    // add hole to the parton list
                      }
                      else
                      {
                          if(outputSource)
                          {
                              // momentum transfer q goes into medium
                              dE = vecq.e();
                              dpx = vecq.px();
                              dpy = vecq.py();
                              dpz = vecq.pz();

                              newSource.type(1);
                              newSource.tau(tau);
                              newSource.x(x);
                              newSource.y(y);
                              newSource.eta(eta);
                              newSource.dE(dE);
                              newSource.dpx(dpx);
                              newSource.dpy(dpy);
                              newSource.dpz(dpz);
                              newSource.vx(vx);
                              newSource.vy(vy);
                              newSource.vz(vz);

                              slist->push_back(newSource);
                          }
                      }
                  }
                  else if( elastic_qg == 1 )   // gluon recoil
                  {
                      // thermal gluon to be recoiled
                      vecThermal = elastic->getThermalMomentum(vecq, T, -1, random);
                      vecRecoilRest = vecq + vecThermal;

                      pxRecoil = vx*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecRecoilRest.px() 
                        + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.py()
                        + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.pz();
                      pyRecoil = vy*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecRecoilRest.py()
                        + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.px()
                        + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.pz();
                      pzRecoil = vz*gamma*vecRecoilRest.pAbs() 
                        + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecRecoilRest.pz()
                        + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.px()
                        + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.py();
                      
                      // momentum four-vector in lab frame
                      vecRecoil.px(pxRecoil);
                      vecRecoil.py(pyRecoil);
                      vecRecoil.pz(pzRecoil);
                      vecRecoil.e(sqrt(pxRecoil*pxRecoil +
                                       pyRecoil*pyRecoil +
                                       pzRecoil*pzRecoil));

                      if( vecRecoilRest.pAbs() > recoilCut )
                      {
                          newOne.p(vecRecoil);
                          newOne.id(21);
                          newOne.mass(0.);
                          newOne.col(700000000+counter);
                          newOne.acol(800000000+counter);

                          newOne.x(x);
                          newOne.y(y);
                          newOne.z(z);
                          newOne.frozen(0);
                          newOne.recoil(1);
                          newOne.init_counts();//Rouz : initialize counts
                          //Sangyong's addition for formation time of radiation
                          newOne.daughter_of(-1); // daughter of nobody
                          newOne.mother_of(-1);  // mother of nobody
                          newOne.p_at_split(newOne.p());

                          plist->push_back(newOne);    // add recoil gluon to the parton list


                          pxThermalLab = vx*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecThermal.px()
                                    + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.py()
                                    + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.pz();
                          pyThermalLab = vy*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecThermal.py()
                                    + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.px()
                                    + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.pz();
                          pzThermalLab = vz*gamma*vecThermal.pAbs() 
                                    + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecThermal.pz()
                                    + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.px()
                                    + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.py();

                          // momentum four-vector in lab frame
                          vecThermalLab.px(pxThermalLab);
                          vecThermalLab.py(pyThermalLab);
                          vecThermalLab.pz(pzThermalLab);
                          vecThermalLab.e(sqrt(pxThermalLab*pxThermalLab +
                                                pyThermalLab*pyThermalLab +
                                                pzThermalLab*pzThermalLab));

                          newHole.p(vecThermalLab);
                          newHole.id(21);
                          newHole.mass(0.);
                          newHole.col(900000000+counter);
                          newHole.acol(1000000000+counter);

                          newHole.x(x);
                          newHole.y(y);
                          newHole.z(z);
                          newHole.frozen(0);
                          newHole.recoil(-1);

                          //Sangyong's addition for formation time of radiation
                          newHole.daughter_of(-1); // daughter of nobody
                          newHole.mother_of(-1);  // mother of nobody
                          newHole.p_at_split(newOne.p());

                          plist->push_back(newHole);    // add hole to the parton list
                      }
                      else
                      {
                          if(outputSource)
                          {
                              // momentum transfer q goes into medium
                              dE = vecq.e();
                              dpx = vecq.px();
                              dpy = vecq.py();
                              dpz = vecq.pz();

                              newSource.type(1);
                              newSource.tau(tau);
                              newSource.x(x);
                              newSource.y(y);
                              newSource.eta(eta);
                              newSource.dE(dE);
                              newSource.dpx(dpx);
                              newSource.dpy(dpy);
                              newSource.dpz(dpz);
                              newSource.vx(vx);
                              newSource.vy(vy);
                              newSource.vz(vz);

                              slist->push_back(newSource);
                          }
                      }
                  }
              }
              else
              {
                  if(outputSource)
                  {
                      // momentum transfer q goes into medium
                      dE = vecq.e();
                      dpx = vecq.px();
                      dpy = vecq.py();
                      dpz = vecq.pz();

                      newSource.type(1);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
              }
          }
      }
      if ( id == 21 ) // if parton is a gluon, let it evolve like a gluon
      {
          // in the following decide what to do in this time step
          // Again, modified for allowing elastic collisions at momenta lower than AMYpCut. -CFY
          // Sangyong's addition:
          // in_coherence is the new flag
          // added if(in_coherence==0) to radiate*
          if (doRadiative == 1 and doElastic == 1 )  //radiative + elastic
              {
                  if(pRest/T>AMYpCut)
                  {
                      totalGluonProb = dt/gamma*(n.Gamma_ggg)+dt/gamma*(n.Gamma_gqq)
                                     + dt/gamma*(gqRate)+dt/gamma*(ggRate)
                                     + dt/gamma*conversiongq; 
                  }
                  else
                  {
                      totalGluonProb = dt/gamma*(gqRate)+dt/gamma*(ggRate)
                   + dt/gamma*conversiongq;
                  }
                  if (totalGluonProb>1)
                  {
                      cout << "MARTINI:WARNING! total probability for gluon to undergo process=" 
                          << totalGluonProb << " > 1. Reduce time step to cure this." << endl; 
                      cout << " pRest: "<< pRest<< " source : "<<plist->at(i).source()<<endl;
                  }
                  if ( random->genrand64_real1() < totalGluonProb ) // check if something happens with the quark
                  {
                      if(pRest/T>AMYpCut)
                      {
                          // Sangyong's addition:
                          // in_coherence is the new flag
                          // added if(in_coherence==0) to radiate* and convert*
                          // this includes coherence effects in radiations, but not in conversions
                          // aways do elastic - elastic was lumped together with the conversions
                          // I am modifying this to use one random number
                          // we have
                          // totalGluonProb = 
                          // dt/gamma*(n.Gamma_ggg)
                          // +dt/gamma*(n.Gamma_gqq)
                          // +dt/gamma*(gqRate)
                          // +dt/gamma*(ggRate)
                          // +dt/gamma*conversiongq; 
                          
                          randx = random->genrand64_real1();
                          
                          AccProb = 0.0;
                          NextProb = (dt/gamma*(gqRate))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                              elastic_gq = 1;

                          AccProb += NextProb;
                          NextProb = (dt/gamma*(ggRate))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                              elastic_gg = 1;

                          AccProb += NextProb;
                          NextProb = dt/gamma*(n.Gamma_ggg)/totalGluonProb;
                          if(AccProb <= randx && randx < (AccProb + NextProb))
                          {
                              if(in_coherence==0) 
                                  radiate_ggg = 1;
                          }
                          
                          AccProb += NextProb;
                          NextProb = dt/gamma*(n.Gamma_gqq)/totalGluonProb;
                          if(AccProb <= randx && randx < (AccProb + NextProb))
                          {
                              if(in_coherence==0) 
                                  radiate_gqq = 1;
                          }
                          
                          AccProb += NextProb;
                          NextProb = (dt/gamma*conversiongq)/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                          {
                              if(in_coherence==0 or plist->at(i).daughter() != -1)
                                  convert_gluon = 1;
                          }
                      }// pRest/T > AMYpCut
                      else
                      {
                      randx = random->genrand64_real1();

                          AccProb = 0.;
                          NextProb = (dt/gamma*(gqRate))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                              elastic_gq = 1;

                          AccProb += NextProb;
                          NextProb = (dt/gamma*(ggRate))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                              elastic_gg = 1;

                          AccProb += NextProb;
                          NextProb = (dt/gamma*(conversiongq))/totalGluonProb; 
                          if( AccProb <= randx && randx < (AccProb + NextProb) )
                          {
                              if(in_coherence==0 or plist->at(i).daughter() != -1)
                                  convert_gluon = 1;
                          }
                      }
                  }
                  else 
                  {
                      continue;
                  }
              }
              else if (doRadiative == 1 and doElastic == 0 and pRest/T>AMYpCut) //radiative
            // no sense implementing coherence here
              {
                  totalGluonProb = dt/gamma*(n.Gamma_ggg)+dt/gamma*(n.Gamma_gqq); 
                  if (totalGluonProb>1) 
                      cout << "MARTINI:WARNING! total probability for gluon to undergo process=" 
                           << totalGluonProb << " > 1. Reduce time step to cure this." << endl; 
                  if ( random->genrand64_real1() < totalGluonProb ) // check if something happens with the quark
                  {
                      if ( random->genrand64_real1() < n.Gamma_ggg/(n.Gamma_ggg+n.Gamma_gqq) ) 
                          radiate_ggg=1; 
                      else 
                          radiate_gqq=1; 
                  }
                  else 
                  {
                      continue;
                  }
              }// rad only
              else if (doRadiative == 0 and doElastic == 1 ) //elastic
              // no sense implementing coherence here
              {
                  totalGluonProb = dt/gamma*(gqRate)+dt/gamma*(ggRate)
                     + dt/gamma*conversiongq; // dt/gamma is dt_rest-frame

                  if (totalGluonProb>1) 
                      cout << "MARTINI:WARNING! total probability for quark to undergo process=" 
                           << totalGluonProb << " > 1. Reduce time step to cure this." << endl; 
                  if ( random->genrand64_real1() < totalGluonProb )  // check if something happens
                  {
                      if ( random->genrand64_real1() < conversiongq/(conversiongq+gqRate+ggRate) )
                      {
                          convert_gluon = 1;
                      }
                      else
                      {
                          if ( random->genrand64_real1() < gqRate/(gqRate+ggRate) )
                          {
                              elastic_gq = 1;
                              elastic_gg = 0;
                          }
                          else
                          {
                              elastic_gg = 1;
                              elastic_gq = 0;
                          }
                      }
                  }
                  else 
                  {
                      continue;
                  }
              }
            else 
                  continue;
          
          // now do what has been decided before
          if(radiate_ggg == 1) // g -> gg
          {
              if(pRest/T>AMYpCut) 
              {
                  f = rates->findValuesRejection
                    (pRest/T, T, alpha_s_rad, Nf, random, import, 3);   // do process 3

                  delp = f.y*T;
                  if(delp > pRest) continue;

                  p = (pRest - delp)*boostBack;
                  newvecp=p/vecp.pAbs()*vecp; // new gluon momentum
                  plist->at(i).p(newvecp);    // change quark's momentum to the new one

                  if(delp > pRest - pCut)
                  {
                      plist->at(i).frozen(2);
                      // exclude soft recoil by assigning id == 0
                      if(plist->at(i).recoil() == -1) plist->at(i).id(0);
                  }
                  if(delp > pCut)
                  {
                      newOne.p(delp*vecp.px()/vecp.pAbs()*boostBack, // emitted gluon's momentum (coll) 
                               delp*vecp.py()/vecp.pAbs()*boostBack,
                               delp*vecp.pz()/vecp.pAbs()*boostBack); 
                      acol = plist->at(i).acol();            // the gluon's anti-color
                      plist->at(i).increment_hard_radiation(); // Rouz: increment hard radiation
                      plist->at(i).acol(1100000000+counter); // set the first new gluon's anti-color to a new one 
                      newOne.id(21);                         // second new particle is a gluon too
                      newOne.col(1100000000+counter);        // set the second gluon's color to the new color 
                      newOne.acol(acol);                     // set second gluon's anti-color to the original one
                      newOne.x(x);
                      newOne.y(y);
                      newOne.z(z);
                      newOne.mass(0.);  // set second gluon's mass to zero
                      newOne.frozen(0);
                      newOne.source(12); //  RY: gluon emitted from a gluon. All photons produced in evolve(...) have itsSource=1
                      newOne.recoil(plist->at(i).recoil());
                      newOne.init_counts();//Rouz:: start the counts for this parton
                      
                      //Sangyong's addition for formation time of radiation
                      daughter = plist->size(); // the position of daughter in the plist when push_back'ed
                      newOne.daughter_of(i); // daughter of the current parton
                      newOne.mother_of(-1);  // mother of nobody
                      newOne.p_at_split(newOne.p()); // momentum at the split point
              
                      plist->at(i).daughter_of(-1); // daughter of nobody
                      plist->at(i).mother_of(daughter);
                      plist->at(i).p_at_split(plist->at(i).p());
                      
                      plist->push_back(newOne);                  // add the second gluon to the list of partons
                  }
                  if ( delp < pCut)
                    plist->at(i).increment_soft_radiation();

                  if(outputSource)
                  {
                      // if delp < 0 or delp > pRest, jet energy gain from thermal medium.
                      // dE is negative but the direction is same as mother jet.
                      // if 0 < delp < pCut or pRest - pCut < delp < pRest, 
                      // radiated jet particle is absorbed into thermal medium.
                      if(delp < pCut)
                      {
                          dE = delp;
                          dpx = fabs(delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(delp)/vecp.pAbs()*vecp.pz();
                      }
                      else if (delp > pRest - pCut)
                      {
                          dE = pRest-delp;
                          dpx = fabs(pRest-delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(pRest-delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(pRest-delp)/vecp.pAbs()*vecp.pz();
                      }
                      newSource.type(2);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
              }
          }
          //pRest = p * gamma * (1.-beta*cosPhi); // boost p to fluid cell's rest frame
          if( radiate_gqq == 1 ) // g -> qq
          {
            if(pRest/T>AMYpCut) 
              {
                  f = rates->findValuesRejection(pRest/T, T, alpha_s_rad, Nf, random, import, 2); // do process 2
                  delp = f.y*T;
                  if(delp > pRest or delp < 0.)
                      continue;

                  p = (pRest - delp)*boostBack;
                  // choose if it is a u-ubar, d-dbar or s-sbar pair:
                  double r = random->genrand64_real1();
                  // note that when using general N_f this has to be changed
                  if ( r < 0.33 ) newID=1;
                  else if  ( r < 0.66 ) newID=2;
                  else newID=3;
                  double mass;
                  if ( newID < 3 ) mass=0.33;  // set the quark's mass
                  else mass = 0.5;         // (pythia needs that for fragmentation) 
                  
                  col = plist->at(i).col();   // gluon's color
                  acol = plist->at(i).acol(); // gluon's anti-color 
                  int col_m, acol_m, col_d, acol_d;
                  double matter = random->genrand64_real1();
                  if (matter<0.5) //flip matter and antimatter
                  {
                      newID *= -1;
                      col_m  = 0;
                      acol_m = acol;
                      col_d  = col;
                      acol_d = 0;
                  }
                  else
                  {
                      col_m  = col;
                      acol_m = 0;
                      col_d  = 0;
                      acol_d = acol;
                  }
                  
                  newvecp = p/vecp.pAbs()*vecp; // new quark momentum
                  plist->at(i).id(newID);     // turn gluon into quarks
                  plist->at(i).col(col_m);    // set quark's color
                  plist->at(i).acol(acol_m);  // set quark's anti-color to zero
                  plist->at(i).mass(mass);    // set quark's mass
                  plist->at(i).p(newvecp);    // change quark's momentum to the new one
                  plist->at(i).source(13);    // this fermion now is the result of gluon splitting -RY.
                  if(p < pCut)
                  {
                       plist->at(i).frozen(2);
                       // exclude soft recoil by assigning id == 0
                       if(plist->at(i).recoil() == -1)
                          plist->at(i).id(0);
                  }
                  plist->at(i).init_counts(); // Rouz: initialize counts for both since 
                  newOne.init_counts();       // this is sort of a conversion process
                  newOne.p(delp*vecp.px()/vecp.pAbs()*boostBack,   // second quark's momentum (collinear)
                         delp*vecp.py()/vecp.pAbs()*boostBack,
                         delp*vecp.pz()/vecp.pAbs()*boostBack); 
                  newOne.id(-newID);                                // second new particle is an anti-quark (minus-sign)
                  newOne.mass(mass);                                // set anti-quark's mass  
                  newOne.col(col_d);                                    // set anti-quark's color to zero
                  newOne.acol(acol_d);                                // set anti-quark's anti-color to gluon's anti-color
                  newOne.x(x);
                  newOne.y(y);
                  newOne.z(z);
                  if(delp > pCut)
                      newOne.frozen(0);
                  else
                      newOne.frozen(2);
                  newOne.source(13); // RY: this fermion is also the result of gluon splitting 
                  newOne.recoil(plist->at(i).recoil());

                  //Sangyong's addition for formation time of radiation
                  daughter = plist->size(); // the position of daughter in the plist when push_back'ed
                  newOne.daughter_of(i);    // daughter of the current parton
                  newOne.mother_of(-1);     // mother of nobody
                  newOne.p_at_split(newOne.p()); // momentum at the split point

                  plist->at(i).daughter_of(-1); // daughter of nobody
                  plist->at(i).mother_of(daughter);
                  plist->at(i).p_at_split(plist->at(i).p());

                  plist->push_back(newOne); // add the second quark to the list of partons
                  if(outputSource)
                  {
                      // if delp < 0 or delp > pRest, jet energy gain from thermal medium.
                      // dE is negative but the direction is same as mother jet.
                      // if 0 < delp < pCut or pRest - pCut < delp < pRest, 
                      // radiated jet particle is absorbed into thermal medium.
                      if(delp < pCut)
                      {
                          dE = delp;
                          dpx = fabs(delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(delp)/vecp.pAbs()*vecp.pz();
                      }
                      else if (delp > pRest - pCut)
                      {
                          dE = pRest-delp;
                          dpx = fabs(pRest-delp)/vecp.pAbs()*vecp.px();
                          dpy = fabs(pRest-delp)/vecp.pAbs()*vecp.py();
                          dpz = fabs(pRest-delp)/vecp.pAbs()*vecp.pz();
                      }
                      newSource.type(2);
                      newSource.tau(tau);
                      newSource.x(x);
                      newSource.y(y);
                      newSource.eta(eta);
                      newSource.dE(dE);
                      newSource.dpx(dpx);
                      newSource.dpy(dpy);
                      newSource.dpz(dpz);
                      newSource.vx(vx);
                      newSource.vy(vy);
                      newSource.vz(vz);

                      slist->push_back(newSource);
                  }
              }
          }
          if ( convert_gluon == 1 ) // do g->q
          {
              if ( random->genrand64_real1() < 0.5 ) 
                GluonToQuark = 1;
              else 
                GluonToQuark = 0;
         
              double r = random->genrand64_real1();
              if (r<0.33) newID=1;
              else if  (r<0.66) newID=2;
              else newID=3;
              
              col = plist->at(i).col();   // color of original gluon  
              acol = plist->at(i).acol(); // anti-color

              if (GluonToQuark == 1)
              {
                plist->at(i).acol(0);     // make the gluon a quark
                newOne.acol(acol);        // let the thermal gluon have the gluon's previous anti-color
                newOne.col(1200000000+counter);
                plist->at(i).id(newID);   // convert to quark
              }
              else
              {
                plist->at(i).col(0); // make the gluon an anti-quark
                newOne.col(col);     // let the thermal gluon have the gluon's previous color
                newOne.acol(1200000000+counter);
                plist->at(i).id(-newID);  // convert to anti-quark
              }
              if (newID == 1 or newID == 2) 
                plist->at(i).mass(0.33);           
              else
                plist->at(i).mass(0.55);
              plist->at(i).source(15);
              plist->at(i).init_counts();// Rouz: initialize counts 
              newOne.init_counts();      // Rouz: initialize counts
              newOne.id(21); // add a thermal gluon
              newOne.mass(0.); 
              newOne.frozen(2);
              newOne.x(x);
              newOne.y(y);
              newOne.z(z);
              newOne.p(random->thermal(T, -1)); // sample momentum from Bose distribution at temperature T
              newOne.source(1);                 // 
              newOne.recoil(plist->at(i).recoil());
              
              //Sangyong's addition for formation time of radiation
              newOne.daughter_of(-1); // daughter of nobody
              newOne.mother_of(-1);  // mother of nobody
              newOne.p_at_split(newOne.p()); // momentum at the conversion point
              
              plist->push_back(newOne); 

              double mass;
              if (newID<3) mass=0.33;  // set the quark's mass
              else mass = 0.5;         // (pythia needs that for fragmentation) 

              if (GluonToQuark == 1)
              {
                newOne.acol(1200000000+counter); // let the thermal anti-quark have the right anti-color
                newOne.col(0);
                newOne.id(-newID);              // add a thermal anti-quark
              }
              else
              {
                newOne.col(1200000000+counter);  // let the thermal quark have the right color
                newOne.acol(0);
                newOne.id(newID);  // add a thermal quark
              }

              newOne.mass(mass); 
              newOne.frozen(0);
              newOne.x(x);
              newOne.y(y);
              newOne.z(z);
              newOne.p(random->thermal(T, 1)); // sample momentum from Fermi distribution at temperature T
              newOne.source(1);               // RY: gluon to fermion conversion. 
              newOne.recoil(plist->at(i).recoil());

              //Sangyong's addition for formation time of radiation
              newOne.daughter_of(-1); // daughter of nobody
              newOne.mother_of(-1);  // mother of nobody
              newOne.p_at_split(newOne.p()); // momentum at the conversion point
              
              plist->push_back(newOne); 

              continue;  //prevent the new gluon from interacting again in this time step 
            }
            else if ( elastic_gq == 1 ) // do gq->gq
            {
              for(;;)
              {
                  omega = elastic->findValuesRejection(pRest/T, T, alpha_s_elas, Nf, random, import, 2); // this is omega/T
                  if(fabs(omega) < pRest/T) break;
              }

              if( transferTransverseMomentum == 1 and omega < 800. )//&& fabs(eta) < 4.)
              {
                  int repeat = 1;
                  double k_min;
                  do
                  {
                      q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s_elas, Nf, random, import, 2);
                      k_min = (q-omega)/2.;
                      if(k_min < 20.) repeat = 0;
                  } while (repeat);
              }
              else
              {
                  q = omega;
              }
            }
            else if ( elastic_gg == 1 ) // do gg->gg
            {
              for(;;)
              {
                  omega = elastic->findValuesRejection(pRest/T, T, alpha_s_elas, Nf, random, import, 4); // this is omega/T
                  if(fabs(omega) < pRest/T) break;
              }

              if( transferTransverseMomentum == 1 && omega < 800. )//&& fabs(eta) < 4.)
                {
                    int repeat = 1;
                    double k_min;
                    do
                    {
                        q = elastic->findValuesRejectionOmegaQ(pRest/T, omega, T, alpha_s_elas, Nf, random, import, 4);
                        k_min = (q-omega)/2.;
                        if(k_min < 20.) repeat = 0;
                    } while (repeat);
                }
                else
                {
                    q = omega;
                }
            }
            if ( elastic_gq == 1 || elastic_gg == 1 ) // in both cases ( gq->gq and gg->gg ) change the momentum 
            {
                newvecpRest=elastic->getNewMomentum(vecpRest, omega*T, q*T, random); // takes everything in GeV

                if (beta > 1e-15) // was beta == 0 but must avoid equality operations with doubles -RY
                {
                    newpx = vx*gamma*newvecpRest.pAbs() 
                      + (1.+(gamma-1.)*vx*vx/(beta*beta))*newvecpRest.px() 
                      + (gamma-1.)*vx*vy/(beta*beta)*newvecpRest.py()
                      + (gamma-1.)*vx*vz/(beta*beta)*newvecpRest.pz();
                    newpy = vy*gamma*newvecpRest.pAbs() 
                      + (1.+(gamma-1.)*vy*vy/(beta*beta))*newvecpRest.py()
                      + (gamma-1.)*vx*vy/(beta*beta)*newvecpRest.px()
                      + (gamma-1.)*vy*vz/(beta*beta)*newvecpRest.pz();
                    newpz = vz*gamma*newvecpRest.pAbs() 
                      + (1.+(gamma-1.)*vz*vz/(beta*beta))*newvecpRest.pz()
                      + (gamma-1.)*vx*vz/(beta*beta)*newvecpRest.px()
                      + (gamma-1.)*vy*vz/(beta*beta)*newvecpRest.py();
                    

                    newvecp.px(newpx);
                    newvecp.py(newpy);
                    newvecp.pz(newpz);
                    newvecp.e(sqrt(newpx*newpx+newpy*newpy+newpz*newpz));
                }
                else
                    newvecp = newvecpRest;
               
                double newpRest = newvecpRest.pAbs();
                if(newpRest < pCut)
                {
                    if(outputSource)
                    {
                        // output energy deposition from jet particle absorption
                        dE = newvecpRest.pAbs();
                        dpx = newvecpRest.px();
                        dpy = newvecpRest.py();
                        dpz = newvecpRest.pz();

                        newSource.type(1);
                        newSource.tau(tau);
                        newSource.x(x);
                        newSource.y(y);
                        newSource.eta(eta);
                        newSource.dE(dE);
                        newSource.dpx(dpx);
                        newSource.dpy(dpy);
                        newSource.dpz(dpz);
                        newSource.vx(vx);
                        newSource.vy(vy);
                        newSource.vz(vz);

                        slist->push_back(newSource);
                    }
                    plist->at(i).frozen(2);
                    //newvecp *= pdummy/newvecp.pAbs();
                    // exclude soft recoil by assigning id == 0
                    if(plist->at(i).recoil() == -1) plist->at(i).id(0);
                }

                plist->at(i).p(newvecp);           // change quark's momentum to the new one

                vecq = vecpRest - newvecpRest;             // momentum transfer q
                if ( vecq.pAbs() < pCut)//Rouz: elastic scattering counts for the gluon
                    plist->at(i).increment_soft_elastic();
                if ( vecq.pAbs() > pCut)
                    plist->at(i).increment_hard_elastic();

                if(getRecoil && omega > 0. && plist->at(i).recoil() != -1)
                {
                   if ( elastic_gq == 1 )  // quark recoil
                   {
                       // thermal quark to be recoiled
                       vecThermal = elastic->getThermalMomentum(vecq, T, 1, random);
                       vecRecoilRest = vecq + vecThermal;

                       pxRecoil = vx*gamma*vecRecoilRest.pAbs() 
                         + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecRecoilRest.px() 
                         + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.py()
                         + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.pz();
                       pyRecoil = vy*gamma*vecRecoilRest.pAbs() 
                         + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecRecoilRest.py()
                         + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.px()
                         + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.pz();
                       pzRecoil = vz*gamma*vecRecoilRest.pAbs() 
                         + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecRecoilRest.pz()
                         + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.px()
                         + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.py();
                       
                       // momentum four-vector in lab frame
                       vecRecoil.px(pxRecoil);
                       vecRecoil.py(pyRecoil);
                       vecRecoil.pz(pzRecoil);
                       vecRecoil.e(sqrt(pxRecoil*pxRecoil +
                                        pyRecoil*pyRecoil +
                                        pzRecoil*pzRecoil));

                       if( vecRecoilRest.pAbs() > recoilCut )
                       {
                           col = plist->at(i).col();            // color of original gluon
                           acol = plist->at(i).acol();          // anti-color of original gluol

                           double r = random->genrand64_real1();
                           // note that when using general N_f this has to be changed
                           if (r<0.33) newID=1;
                           else if (r<0.66) newID=2;
                           else newID=3;
                           double mass;
                           if (newID<3) mass=0.33;    // set the quark's mass
                           else mass = 0.5;           // (pythia needs that for fragmentation) 

                            newOne.p(vecRecoil);

                            newOne.mass(mass);
                            double matter = random->genrand64_real1();

                            if(matter < 0.5)   // gq -> gq
                            {
                                newOne.id(newID);
                                newOne.col(1300000000+counter);
                                newOne.acol(0);

                                newHole.id(newID);
                                newHole.col(1400000000+counter);
                                newHole.acol(0);
                            }
                        else  // g qbar -> g qbar
                        {
                                newOne.id(-newID);
                                newOne.col(0);
                                newOne.acol(1300000000+counter);

                                newHole.id(-newID);
                                newHole.col(0);
                                newHole.acol(1400000000+counter);
                        }

                            newOne.x(x);
                            newOne.y(y);
                            newOne.z(z);
                            newOne.frozen(0);
                            newOne.recoil(1);
                            newOne.init_counts(); // Rouz: initialize counts for the recoil particle off gluon 
                            //Sangyong's addition for formation time of radiation
                            newOne.daughter_of(-1); // daughter of nobody
                            newOne.mother_of(-1);  // mother of nobody
                            newOne.p_at_split(newOne.p()); // momentum at the conversion point

                            plist->push_back(newOne);  // add the gluon to the list of partons if (f.y>AMYpCut)

                            pxThermalLab = vx*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecThermal.px()
                                      + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.py()
                                      + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.pz();
                            pyThermalLab = vy*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecThermal.py()
                                      + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.px()
                                      + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.pz();
                            pzThermalLab = vz*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecThermal.pz()
                                      + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.px()
                                      + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.py();


                            // momentum four-vector in lab frame
                            vecThermalLab.px(pxThermalLab);
                            vecThermalLab.py(pyThermalLab);
                            vecThermalLab.pz(pzThermalLab);
                            vecThermalLab.e(sqrt(pxThermalLab*pxThermalLab +
                                                  pyThermalLab*pyThermalLab +
                                                  pzThermalLab*pzThermalLab));

                            newHole.p(vecThermalLab);
                            newHole.mass(mass);

                            newHole.x(x);
                            newHole.y(y);
                            newHole.z(z);
                            newHole.frozen(0);
                            newHole.recoil(-1);
                            
                            //Sangyong's addition for formation time of radiation
                            newHole.daughter_of(-1); // daughter of nobody
                            newHole.mother_of(-1);  // mother of nobody
                            newHole.p_at_split(newOne.p()); // momentum at the conversion point

                            plist->push_back(newHole);    // add hole to the parton list
                       }
                       else
                       {
                           if(outputSource)
                           {
                               // momentum transfer q goes into medium
                               dE = vecq.e();
                               dpx = vecq.px();
                               dpy = vecq.py();
                               dpz = vecq.pz();

                               newSource.type(1);
                               newSource.tau(tau);
                               newSource.x(x);
                               newSource.y(y);
                               newSource.eta(eta);
                               newSource.dE(dE);
                               newSource.dpx(dpx);
                               newSource.dpy(dpy);
                               newSource.dpz(dpz);
                               newSource.vx(vx);
                               newSource.vy(vy);
                               newSource.vz(vz);

                               slist->push_back(newSource);
                           }
                       }
                   }
                   else if ( elastic_gg == 1 )  // gluon recoil
                   {
                        // thermal gluon to be recoiled
                        vecThermal = elastic->getThermalMomentum(vecq, T, -1, random);
                        vecRecoilRest = vecq + vecThermal;

                        pxRecoil = vx*gamma*vecRecoilRest.pAbs() 
                          + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecRecoilRest.px() 
                          + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.py()
                          + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.pz();
                        pyRecoil = vy*gamma*vecRecoilRest.pAbs() 
                          + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecRecoilRest.py()
                          + (gamma-1.)*vx*vy/(beta*beta)*vecRecoilRest.px()
                          + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.pz();
                        pzRecoil = vz*gamma*vecRecoilRest.pAbs() 
                          + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecRecoilRest.pz()
                          + (gamma-1.)*vx*vz/(beta*beta)*vecRecoilRest.px()
                          + (gamma-1.)*vy*vz/(beta*beta)*vecRecoilRest.py();
                        
                        // momentum four-vector in lab frame
                        vecRecoil.px(pxRecoil);
                        vecRecoil.py(pyRecoil);
                        vecRecoil.pz(pzRecoil);
                        vecRecoil.e(sqrt(pxRecoil*pxRecoil +
                                         pyRecoil*pyRecoil +
                                         pzRecoil*pzRecoil));

                      if( vecRecoilRest.pAbs() > recoilCut )
                      {
                            newOne.p(vecRecoil);
                            newOne.id(21);
                            newOne.mass(0.);
                            newOne.col(1500000000+counter);
                            newOne.acol(1600000000+counter);
                            newOne.x(x);
                            newOne.y(y);
                            newOne.z(z);
                            newOne.frozen(0);
                            newOne.recoil(1);
                            newOne.init_counts(); // Rouz: recoil initialize counts
                            //Sangyong's addition for formation time of radiation
                            newOne.daughter_of(-1); // daughter of nobody
                            newOne.mother_of(-1);  // mother of nobody
                            newOne.p_at_split(newOne.p()); // momentum at the conversion point

                            plist->push_back(newOne);    // add recoil parton to the parton list

                            pxThermalLab = vx*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vx*vx/(beta*beta))*vecThermal.px()
                                      + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.py()
                                      + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.pz();
                            pyThermalLab = vy*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vy*vy/(beta*beta))*vecThermal.py()
                                      + (gamma-1.)*vx*vy/(beta*beta)*vecThermal.px()
                                      + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.pz();
                            pzThermalLab = vz*gamma*vecThermal.pAbs() 
                                      + (1.+(gamma-1.)*vz*vz/(beta*beta))*vecThermal.pz()
                                      + (gamma-1.)*vx*vz/(beta*beta)*vecThermal.px()
                                      + (gamma-1.)*vy*vz/(beta*beta)*vecThermal.py();


                            // momentum four-vector in lab frame
                            vecThermalLab.px(pxThermalLab);
                            vecThermalLab.py(pyThermalLab);
                            vecThermalLab.pz(pzThermalLab);
                            vecThermalLab.e(sqrt(pxThermalLab*pxThermalLab +
                                                  pyThermalLab*pyThermalLab +
                                                  pzThermalLab*pzThermalLab));

                            newHole.p(vecThermalLab);

                            newHole.id(21);
                            newHole.mass(0.);
                            newOne.col(1700000000+counter);
                            newOne.acol(1800000000+counter);

                            newHole.x(x);
                            newHole.y(y);
                            newHole.z(z);
                            newHole.frozen(0);
                            newHole.recoil(-1);

                            //Sangyong's addition for formation time of radiation
                            newHole.daughter_of(-1); // daughter of nobody
                            newHole.mother_of(-1);  // mother of nobody
                            newHole.p_at_split(newOne.p()); // momentum at the conversion point

                            plist->push_back(newHole);    // add hole to the parton list
                        }
                        else
                        {
                            if(outputSource)
                            {
                                // momentum transfer q goes into medium
                                dE = vecq.e();
                                dpx = vecq.px();
                                dpy = vecq.py();
                                dpz = vecq.pz();

                                newSource.type(1);
                                newSource.tau(tau);
                                newSource.x(x);
                                newSource.y(y);
                                newSource.eta(eta);
                                newSource.dE(dE);
                                newSource.dpx(dpx);
                                newSource.dpy(dpy);
                                newSource.dpz(dpz);
                                newSource.vx(vx);
                                newSource.vy(vy);
                                newSource.vz(vz);

                                slist->push_back(newSource);
                            }
                        }
                    }
                }
                else
                {
                    if(outputSource)
                    {
                        // momentum transfer q goes into medium
                        dE = vecq.e();
                        dpx = vecq.px();
                        dpy = vecq.py();
                        dpz = vecq.pz();

                        newSource.type(1);
                        newSource.tau(tau);
                        newSource.x(x);
                        newSource.y(y);
                        newSource.eta(eta);
                        newSource.dE(dE);
                        newSource.dpx(dpx);
                        newSource.dpy(dpy);
                        newSource.dpz(dpz);
                        newSource.vx(vx);
                        newSource.vy(vy);
                        newSource.vz(vz);

                        slist->push_back(newSource);
                    }
                }
            }
        }
    } // loop over partons i

    return counter;
}

void MARTINI::sampleTA()
{
    int A = static_cast<int>(glauber->get_Projetile_A());
    //cout << "A=" << A << endl;

    ReturnValue rv, rv2;
    
    for (int i = 0; i < A; i++) // get all nucleon coordinates
    {
        rv = glauber->SampleTARejection(random);
        rv2 = glauber->SampleTARejection(random);
        nucleusA.push_back(rv);
        nucleusB.push_back(rv2);
    }
}  

int MARTINI::generateEvent(vector<Parton> *plist)
{
   if(fullEvent==1)
   {
      // sample Nu (which gives the number of nucleons at radial position r) for both nuclei
      // then lay them on top of each other and sample the number of jet events in a tube with area=inelasticXSec
      
      //cout << " number of binary collisions = " << glauber->TAB() << endl;
      //cout << " A=" << glauber->nucleusA() << endl;
      Parton parton;                               // new Parton object
      double b = glauberImpactParam;
      int posx;
      int posy;
      int n1=0;
      int n2=0;
      double r = 1.2*pow(glauber->get_Projetile_A(), 1./3.);

      nucleusA.clear();
      nucleusB.clear();
      sampleTA();                                  // populate the lists nucleusA and nucleusB with position data of the nucleons
      
      double cellLength = sqrt(inelasticXSec*0.1); // compute cell length in fm (1 mb = 0.1 fm^2)
      double latticeLength = 4.*r;                 // spread the lattice 2r in both (+/-) directions (total length = 4r)
      int ixmax = ceil(latticeLength/cellLength);
      ixmax*=2;                                    // number of cells in the x (and y) direction 
      double xmin = -latticeLength;                // minimal x (and y) on the lattice
      
      int nucALat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus A
      int nucBLat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus B
      int collLat[ixmax][ixmax];                   // lattice with the positions of the interactions -> for information only
      
      int countCollisionsA[100];
      int countCollisionsB[100];

      for (int i=0; i<100; i++)
      {
        countCollisionsA[i]=0;
        countCollisionsB[i]=0;
      }

      for (int i = 0; i < ixmax; i++)              // initialize cells: zero nucleons in every cell
        for (int j = 0; j < ixmax; j++)
        {
           nucALat[i][j] = 0;
           nucBLat[i][j] = 0;
           collLat[i][j] = 0;
        }

      
      //cout <<"size=" << nucleusA.size() << endl;
      for (int i = 0; i < nucleusA.size(); i++)    // fill nucleons into cells
      {
        posx = floor((nucleusA.at(i).x-xmin-b/2.)/cellLength);
        posy = floor((nucleusA.at(i).y-xmin)/cellLength);
        nucALat[posx][posy] += 1;
        
        posx = floor((nucleusB.at(i).x-xmin+b/2.)/cellLength);
        posy = floor((nucleusB.at(i).y-xmin)/cellLength);
        nucBLat[posx][posy] += 1;
      }
      
      int done;
      int numberOfpp;
      int eventNumber=0;
      numberOfpp=0;
      //cout << "jetXSec=" << jetXSec << endl;
      //cout << "inelasticXSec=" << inelasticXSec << endl;
      int Projectile_A=glauber->get_Projetile_A();
      int Target_A=glauber->get_Projetile_A();
      int Projectile_Z = glauber->get_Projetile_Z();
      int Target_Z = glauber->get_Projetile_Z();
      //cout << " probability to generate jet event =" << jetXSec/inelasticXSec << endl;
      for (int ix = 0; ix < ixmax; ix++) 
      {
        for (int iy = 0; iy < ixmax; iy++)
        {
          for (int i = 0; i < nucALat[ix][iy]; i++)
          {
            for (int j = 0; j < nucBLat[ix][iy]; j++)
            {
              //cout << "ix+b/2=" << static_cast<int>(ix+b/2.) << " # there=" << nucALat[static_cast<int>(ix+b/2.)][iy] << endl;
              if ( random->genrand64_real1() < jetXSec/inelasticXSec )
              {
                countCollisionsA[i]+=1;
                countCollisionsB[j]+=1;
                double xPositionInCell = (random->genrand64_real1())*cellLength;
                double yPositionInCell = (random->genrand64_real1())*cellLength;
                // see if colliding nucleons are neutrons
                if(random->genrand64_real1() > (double)Projectile_Z/(double)Projectile_A) 
                  n1=2112;
                else 
                  n1=2212;
                if(random->genrand64_real1() > (double)Target_Z/(double)Target_A) 
                  n2=2112;
                else 
                  n2=2212;
                pythia.next();
//                pythia.next(n1,n2);                                                // generate event with pythia
                done = 0;
                //pythia.event.list();
                for (int ip = 0; ip < pythia.event.size(); ++ip) 
                {
                  // if the parton is final, i.e., present after the showers, then put it in the list
                  // currently heavy quarks are taken but will not lose energy! 
                  // only g and u,d,s lose energy while the 3 quarks are assumed to be massless
                  if (pythia.event[ip].isFinal()
                      // && (pythia.event[i].id()<4 || pythia.event[i].id()==21)
                     )
                  {
                    parton.id(pythia.event[ip].id());                     // set parton id
                    parton.status(pythia.event[ip].status());             // set parton status
                    parton.mass(pythia.event[ip].m());                    // set mass
                    if (fixedTemperature==0)
                    {
                      //parton.x(xmin+(ix)*cellLength+cellLength/2.);       // set position
                      //parton.y(xmin+(iy)*cellLength+cellLength/2.);
                      parton.x(xmin+(ix)*cellLength+xPositionInCell);       // set position
                      parton.y(xmin+(iy)*cellLength+yPositionInCell);
                      //parton.xini(xmin+(ix)*cellLength+xPositionInCell);       // set position
                      //parton.yini(xmin+(iy)*cellLength+yPositionInCell);
                      if ( done == 0 )
                      {
                        //cout << "ix=" << ix << " iy=" << iy << " nucli=" << i 
                        //     << " nuclj=" << j << " ip=" << ip << " x=" << parton.x()<< " y=" << parton.y() << endl;
                        posx = floor((parton.x()-xmin)/cellLength);
                        posy = floor((parton.y()-xmin)/cellLength);
                        collLat[posx][posy] += 1;
                        
                        //cout << numberOfpp << endl;
                        numberOfpp++;
                        done = 1;
                      }
                    }
                    else
                    {
                      parton.x(0.);     
                      parton.y(0.);
                    }
                    parton.z(0.);
                    parton.col(pythia.event[ip].col());                     // set color 
                    parton.acol(pythia.event[ip].acol());                 // set anti-color
                    parton.frozen(0);                                     // parton is not frozen (will evolve)
                    parton.p(pythia.event[ip].p());                       // set momentum
                    parton.source(0);                                     // All initial partons have itsSource=0
                  
                    // Sangyong's addition
                    parton.mother_of(-1); // mother of nobody
                    parton.daughter_of(-1); // daughter of nobody
                    parton.p_at_split(parton.p()); // original momentum, don't need it but might as well populate it
                    // Sangyong's addition

                    plist->push_back(parton);                             // add the parton to the main list
                  }
                } 
                eventNumber++;
              }
            }
          }
        }
      }

      //for (int i=0; i<100; i++)
      //{
      //  cout << "collisions for parton i=" << i << " =" << countCollisionsA[i] << endl;
      //  cout << "collisions for parton j=" << i << " =" << countCollisionsB[i] << endl;
      //}

      //cout << "eventNumber = " << eventNumber << endl;
      //cout << "total number of NN =" << numberOfpp << endl;
      //cout << "eventNumber =" << eventNumber << endl;
      totalNNs = eventNumber;
      //output for plot
      ofstream fout1("./output/density.dat",ios::out); 
      ofstream fout2("./output/densityB.dat",ios::out); 
      ofstream fout3("./output/colls.dat",ios::out); 
      ofstream fout4("./output/NN.dat",ios::app); 
      
      for (int i = 0; i < ixmax; i++)
      {
        for (int j = 0; j < ixmax; j++)
        {
          fout1 <<  i << " " << j << " " << nucALat[i][j] << endl;
          if ( j == ixmax-1 ) fout1 << endl;
        }
      }
      
      for (int i = 0; i < ixmax; i++)
      {
        for (int j = 0; j < ixmax; j++)
        {
          fout2 <<  i << " " << j << " " << nucBLat[i][j] << endl;
          if ( j == ixmax-1 ) fout2 << endl;
        }
      }
      
      for (int i = 0; i < ixmax; i++)
      {
        for (int j = 0; j < ixmax; j++)
        {
          fout3 <<  i << " " << j << " " << collLat[i][j] << endl;
          if ( j == ixmax-1 ) fout3 << endl;
        }
      }
      
      fout4 << numberOfpp << endl;
      
      fout1.close();
      fout2.close();
      fout3.close();
      fout4.close();
      //cout << " done generating initial hard collisions" << endl;
      return eventNumber;
    }
    else if (fullEvent==0) // sample only one hard collision per heavy-ion event
    {
        Parton parton;   // new Parton object
        int n1;
        int n2;
        int Projectile_A = glauber->get_Projetile_A();
        int Target_A = glauber->get_Target_A();
        int Projectile_Z = glauber->get_Projetile_Z();
        int Target_Z = glauber->get_Target_Z();
        if(random->genrand64_real1() > (double)Projectile_Z/(double)Projectile_A)
            n1=2112;
        else 
            n1=2212;
        if(random->genrand64_real1() > (double)Target_Z/(double)Target_A)
            n2=2112;
        else 
            n2=2212;
        pythia.next();
//        pythia.next(n1,n2);    // generate event with pythia
        //pythia.event.list(); 
        
        //ReturnValue rv;
        //if (evolution && fixedTemperature==0) 
        if (!allFromCenter)
        {
            if(nbinFromFile == 1)
            {
                //Sample a random element of the list of number of collisions. -CFY 11/2/2010
                double randomr = random->genrand64_real2();
                int randomint = (int)(randomr*((double)Nbin));
                rv.x = binary[randomint][0];
                rv.y = binary[randomint][1];
                cout << "rv.x = " << rv.x << ", rv.y = " << rv.y << endl;
            }
            else if(nbinFromFile == 2 or nbinFromFile == 3 or nbinFromFile == 4)
            {
               double x_sample, y_sample;
               binary_info_ptr->get_sample(&x_sample, &y_sample);
               rv.x = x_sample;
               rv.y = y_sample;
            }
            else if (glauberEnvelope)
                rv=glauber->SamplePAB(random);          // determine position in x-y plane from Glauber model with Metropolis
            else 
                rv=glauber->SamplePABRejection(random); // determine position in x-y plane from Glauber model with rejection method
            //cout << rv.x << " " << rv.y << endl;
        }
        else
        {
            rv.x=0.;
            rv.y=0.;
        }

        for (int i = 0; i < pythia.event.size(); ++i) 
        {
            // if the parton is final, i.e., present after the showers, then put it in the list
            // currently heavy quarks are taken but will not lose energy! 
            // only g and u,d,s lose energy and the 3 quarks are assumed to be massless
            if (pythia.event[i].isFinal()
                // && (pythia.event[i].id()<4 || pythia.event[i].id()==21)
               )
            {
                parton.id(pythia.event[i].id());         // set parton id
                parton.status(pythia.event[i].status()); // set parton status
                parton.mass(pythia.event[i].m());        // set mass
                
                parton.x(rv.x);                          // set position
                parton.y(rv.y);
                parton.z(0.);
                parton.col(pythia.event[i].col());       // set color 
                parton.acol(pythia.event[i].acol());     // set anti-color
                parton.frozen(0);                        // parton is not frozen (will evolve)
                
                parton.p(pythia.event[i].p());           // set momentum
                parton.source(0);                        // All intial partons have itsSource=0
                parton.recoil(0);                        // Hard partons have itsRecoil=0
                        
                // Sangyong's addition
                parton.mother_of(-1);                   // mother of nobody
                parton.daughter_of(-1);                 // daughter of nobody
                parton.p_at_split(parton.p());          // original momentum, don't need it but might as well populate it
                // Sangyong's addition
                
                plist->push_back(parton);               // add the parton to the main list
                // check jettiness of di-quarks:
                // if (abs(parton.id())>2000) 
                //   cout << "id=" << parton.id() << " px=" << parton.p().px() 
                //        << " py=" << parton.p().py() << " pz=" <<  parton.p().pz() << endl;
            }
        } 
      return 1;
    }
}


int MARTINI::fragmentation(vector<Parton> * plist, int recoil, int currentEvent )
{
    if (fragmentationMethod == 1)
    {
        if (fullEvent == 0 )
        {
            // clear the event record to fill in the partons resulting from evolution below
            pythia.event.reset(); 
            int id, col, acol, status;
            int numEvents = 0;
            double mass;
            Vec4 pvec; 
            vector<Parton> * flist = new vector<Parton>;
            
            for(int p=0; p<plist->size(); p++)
                if( plist->at(p).recoil() == recoil and plist->at(p).id() != 0 )
                    flist->push_back(plist->at(p));
            // reconnect color for recoils/holes
            if(recoil)
            {

                //vector<Parton> *templist;
                //templist = new vector<Parton>;
                //templist->insert(templist->begin(), flist->begin(), flist->begin()+10);
                ////templist->insert(templist->begin(), flist->begin(), flist->end());
                //flist->clear();
                //flist = templist;
                
                //Vec4 pvec;
                //int col, acol;
                //int f_imax = flist->size();
                //for ( int f_i=0; f_i<f_imax; f_i++)
                //{
                //    pvec = flist->at(f_i).p();
                //    col = flist->at(f_i).col();
                //    acol = flist->at(f_i).acol();
                //    cout << "hadronization::col/acol = " << col << " " << acol << " p = " << pvec;
                //}
                
                
                //flist->at(0).id(1);
                //flist->at(1).id(21);
                //flist->at(2).id(21);
                //flist->at(3).id(1);
                //flist->at(4).id(21);
                //flist->at(5).id(21);


                cout << setprecision(3);
                // loop over all partons in the sublist and reset color links
                for (int k=0; k<flist->size(); k++)  
                {
                    flist->at(k).col(0);
                    flist->at(k).acol(0);
                }
                
                unsigned it = 0;
                while(it < flist->size())
                {
                    double T = hydroTfinal;
                    Parton newOne;
                    
                    double r, rn;
                    int newID;
                    double mass;

                    bool quarkLoop = true;

                    // search upcoming q/qbar
                    if (flist->at(it).id() > 0 and flist->at(it).id() < 4)
                    {
                        quarkLoop = true;
                    }
                    else if (flist->at(it).id() < 0)
                    {
                        quarkLoop = false;
                    }
                    // if next parton is g, search upcoming q/qbar
                    else if (flist->at(it).id() == 21)
                    {  
                        //cout << "loop begins with gluon" << endl;
                        int pos = -1;
                        unsigned ii = it+1;
                        while(ii < flist->size())
                        {
                            if(flist->at(ii).id() != 21)
                            {
                                pos = ii;
                                if(flist->at(ii).id() > 0 and flist->at(ii).id() < 4)
                                    quarkLoop = true;
                                else
                                    quarkLoop = false;

                                break;
                            }
                            ii++;
                            if(ii == flist->size()) 
                            {
                                //cout << "no more quark to start loop! Add one at the front." << endl;
                                //cout << "[test]::it = " << it << endl;

                                r = random->genrand64_real1();
                                // note that when using general N_f this has to be changed
                                if (r<0.33) newID=1;
                                else if (r<0.66) newID=2;
                                else newID=3;
                                rn=random->genrand64_real1();
                                if(rn < 0.5) newID *= -1;
                                if (fabs(newID<3)) mass=0.33;
                                else mass = 0.5;

                                newOne.id(newID);
                                newOne.mass(mass);
                                newOne.col(0);
                                newOne.acol(0);

                                newOne.p(random->thermal(T, 1));

                                //cout << "put " << newOne.id() << " at " << it+1 << endl;
                                flist->insert(flist->begin()+it, newOne);
                                //for(unsigned i = 0; i != flist->size(); i++) {
                                //   cout << flist->at(i).id() << " ";
                                //}
                                //cout << endl;
                                break;
                            }
                        }
                        //cout << "pos = " << pos << endl;
                        if(pos > 0)
                        {
                            // rotate current gluon and q/qbar in pos
                            //cout << "[1]rotate::from " << pos+1 << " to " << it+1 << endl;
                            rotate( flist->begin()+pos, flist->begin()+pos+1, flist->begin()+it+1 );
                            //for(unsigned i = 0; i != flist->size(); i++) {
                            //cout << flist->at(i).id() << " ";
                            //}
                            //cout << endl;  
                        }
                    }
                    //cout << "[Interm]::Loop is always started with q/qbar " << flist->at(it).id() << " at " << it+1 << endl;
                    //for(unsigned i = 0; i != flist->size(); i++) {
                    //cout << flist->at(i).id() << " ";
                    //}
                    //cout << endl;
                    // Now go to the position where there is non gluon
                    while(it < flist->size())
                    {
                        it++;
                        // If this q/qbar is the last, add qbar/q to close the loop
                        if(it == flist->size()) 
                        {
                            //cout << "gluon ending or this q/qbar is end! Add one" << endl;
                            r = random->genrand64_real1();
                            // note that when using general N_f this has to be changed
                            if (r<0.33) newID=1;
                            else if (r<0.66) newID=2;
                            else newID=3;
                            if(quarkLoop) newID *= -1;
                            if (fabs(newID<3)) mass=0.33;
                            else mass = 0.5;

                            newOne.id(newID);
                            newOne.mass(mass);
                            newOne.col(0);
                            newOne.acol(0);

                            newOne.p(random->thermal(T, 1));

                            //cout << "put " << newOne.id() << " at " << it+1 << endl;
                            flist->insert(flist->begin()+it, newOne);
                            //for(unsigned i = 0; i != flist->size(); i++) {
                            //   cout << flist->at(i).id() << " ";
                            //}
                            //cout << endl;
                            break;
                        }
                        // Position where there is non-gluon
                        if(flist->at(it).id() != 21)
                            break;

                    }
                    //cout << "[Interm]::next q/bar position : " << it+1 << endl;

                    // 1: close the loop if next one is its counterpart (do nothing)
                    // 2: replace next q/qbar with upcoming counter part with next one
                    // 3: create its counter part if none is found
                    unsigned ii = it;  // ii : position of q/qbar for closing the loop
                    while(ii < flist->size())
                    {
                        // If we find q/qbar to close the loop.
                        if( (quarkLoop && flist->at(ii).id() < 0) ||
                            (!quarkLoop && flist->at(ii).id() > 0 && flist->at(ii).id() < 4) )
                        {
                            //cout << "[2]rotate::from " << ii+1 << " to " << it+1 << endl;
                            rotate( flist->begin()+ii, flist->begin()+ii+1, flist->begin()+it+1 );
                            //for(unsigned i = 0; i != flist->size(); i++) {
                            //    cout << flist->at(i).id() << " ";
                            //}
                            //cout << endl;
                            break;
                        }
                        ii++;
                        // If we don't find any q/qbar to close the loop.
                        if(ii == flist->size()) 
                        {
                            //cout << "no more q/qbar to close loop! Add one" << endl;
                            // add q/qbar if nothing is found
                            unsigned iq = it;  // iq : position of q/qbar to add
                            while(iq < flist->size())
                            {
                                if(flist->at(iq).id() != 21 || iq == flist->size())
                                    break;
                                iq++;
                            }

                            r = random->genrand64_real1();
                            // note that when using general N_f this has to be changed
                            if (r<0.33) newID=1;
                            else if (r<0.66) newID=2;
                            else newID=3;
                            if(quarkLoop) newID *= -1;
                            if (fabs(newID<3)) mass=0.33;
                            else mass = 0.5;

                            newOne.id(newID);
                            newOne.mass(mass);
                            newOne.col(0);
                            newOne.acol(0);

                            newOne.p(random->thermal(T, 1));

                            //cout << "put " << newOne.id() << " at " << it+1 << endl;
                            flist->insert(flist->begin()+iq, newOne);

                            //for(unsigned i = 0; i != flist->size(); i++) {
                            //    cout << flist->at(i).id() << " ";
                            //}
                            //cout << endl;
                            break;
                        }
                    }

                    it++;
                    if(it == flist->size()) 
                    {
                        //cout << "process finished!" << endl;
                        break;
                    }
                    
                    //cout << "[finish]::it = " << it+1 << " new flist->size() = " << flist->size() << endl;
                    //if(it > 5000) exit(3);
                }
                // Now assign col/acol accordingly.
                int stringCounter = 100;
                unsigned ip = 0;
                while(ip < flist->size())
                {
                
                    bool quarkLoop = true;
                    if(flist->at(ip).id() < 0) quarkLoop = false;
                    //cout << "quarkLoop = " << quarkLoop << endl;
                
                    if(quarkLoop)
                        flist->at(ip).acol(0);
                    else
                        flist->at(ip).col(0);
                
                    // closeLoop = true : ending a loop
                    bool closeLoop = false;
                    unsigned ii = ip+1;
                    while(ii < flist->size())
                    {
                        double px1, py1, pz1, p1;
                        double px2, py2, pz2, p2;
                        double cosP1P2;
                
                        px1 = flist->at(ii-1).p().px();
                        py1 = flist->at(ii-1).p().py();
                        pz1 = flist->at(ii-1).p().pz();
                        p1 = flist->at(ii-1).p().pAbs();
                
                        px2 = flist->at(ii).p().px();
                        py2 = flist->at(ii).p().py();
                        pz2 = flist->at(ii).p().pz();
                        p2 = flist->at(ii).p().pAbs();
                        cosP1P2 = (px1*px2 + py1*py2 + pz1*pz2)/(p1*p2);
                
                        if(quarkLoop)
                        {
                            flist->at(ii-1).col(stringCounter);
                            flist->at(ii).acol(stringCounter);
                            if(flist->at(ii).id() != 21) closeLoop = true;
                        }
                        else
                        {
                            flist->at(ii).col(stringCounter);
                            flist->at(ii-1).acol(stringCounter);
                            if(flist->at(ii).id() != 21) closeLoop = true;
                        }
                
                        stringCounter++;
                        if(closeLoop)
                        {
                            ip = ii;
                            break;
                        }
                        ii++;
                    }
                    if(quarkLoop)
                        flist->at(ip).col(0);
                    else
                        flist->at(ip).acol(0);
                
                    ip++;
                }


            }
            // end of color reconnection

            int f_imax = flist->size();

            // add all the higher momentum partons:
            // loop over all partons in the main list
            for ( int f_i=0; f_i<f_imax; f_i++)
            {
                id = flist->at(f_i).id();

                status = flist->at(f_i).status();
                if ( status<1 or status>100)   // if status was not set, cure this here
                      status = 62;
                

                pvec = flist->at(f_i).p();
                col = flist->at(f_i).col();
                acol = flist->at(f_i).acol();
                mass = flist->at(f_i).mass();

                pythia.event.append(id,status,col,acol,pvec,mass);
            }
            // does the fragmentation and the rest (if success add 1 to total events)
            if(pythia.forceHadronLevel()) 
                numEvents+=1; 

            delete flist;
            return numEvents;
        }
        else if (fullEvent == 1)
        {
            //pythia.event.list(); 
            pythia.event.reset(); // clear the event record to fill in the partons resulting from evolution below
            int j = 0;
            int id, col, acol, sid, scol, sacol, status, eventNumber;
            int p_imax = plist->size();
            int p_jmax = plist->size();
            int numEvents = 0;
            //int check = 0;
            double mass, smass;
            Vec4 pvec, spvec; 
            // add all the higher momentum partons:
            for ( int p_i=0; p_i<p_imax; p_i++) // loop over all partons in the main list (those above momentum threshold)
            {
                id = plist->at(p_i).id();
                if (id == 22)   // do not feed photons into fragmentation 
                    continue;   // in this way only fragmentation and decay photons appare at the end
                                // photons from initial state and jet-medium interaction will be collected 
                                // before fragmenetation in the main.cpp

                status = plist->at(p_i).status();
                //if (status<1 || status>100) status = 1; // if status was not set, cure this here
                if (status<1 || status>100)   // if status was not set, cure this here
                      status = 62;

                pvec = plist->at(p_i).p();
                col = plist->at(p_i).col();
                acol = plist->at(p_i).acol();
                mass = plist->at(p_i).mass();
                // append( id, status, col, acol, p, mass )
                if (eventNumber==currentEvent)
                {
                    //if(id == 22)
                    //{
                    //  int source = plist->at(p_i).source();
                    //  cout << scientific << setw(16) << setprecision(6) 
                    //       << "photon: status = " << status << " source = " << source
                    //       << ", px = " << pvec.px() << ", py = " << pvec.py() << ", pz = " << pvec.pz() 
                    //       << endl;
                    //}
                    pythia.event.append(id,status,col,acol,pvec,mass); // copy existing parton into event record
                }
            }
            if(pythia.forceHadronLevel()) numEvents+=1; // does the fragmentation and the rest (if success add 1 to total events)
            //if(check>0) pythia.event.list(); 
            //if(check>0) pythia.event.list(); 
            return numEvents;
        }
    } 
    else if (fragmentationMethod == 2)
    {
        pythia.event.reset(); // clear the event record to fill in the partons resulting from evolution below
        int j = 0;
        int i, is, k, l, ix, iy, iz, ixmax, izmax;
        int Ncells;
        int id, col, acol, sid, scol, sacol, status;
        int p_imax = plist->size();
        int p_jmax = plist->size();
        int numEvents = 0;
        //int check = 0;
        HydroInfo hydroInfo;
        double mass, smass, T, x, y, z, t, V;
        Vec4 pvec, spvec; 
        double Ngluons, Nquarks;
        int Igluons, Iquarks;
        int stringCounter, rni;
        double DX, DZ;
        double rn;
        Parton newOne, temp;
        int countStrange = 0;
        int countAStrange = 0;
        int countq = 0;
        int countaq = 0;
        int countQuarkJets;
        int countAntiQuarkJets;
        
        double eps = 1e-15;
        hydroZmax = 40;
        //cout << "doing fragmentation, method 2" << endl;
        if (hydroWhichHydro != 3 && hydroWhichHydro != 4 && hydroWhichHydro != 7 && hydroWhichHydro !=8) 
        {    
            hydroZmax=hydroTauMax;
            hydroDz=hydroDx;
        }
      
        DX = 1.; //fm
        DZ = 1.; //fm
        // these values are for internal consumption of this function
        
        ixmax = floor((2.*hydroXmax)/DX+0.0001);
        izmax = floor((2.*hydroZmax)/DZ+0.0001);//z is really eta

        Ncells = ixmax*ixmax*izmax;

        vector<Parton> ** psublist;               // pointer to array of vector<Parton> objects
        vector<double> ** cvec;

        psublist = new vector<Parton> *[Ncells];  // pointer to array of vector<Parton> objects
        cvec = new vector<double> *[Ncells]; 

        //cout << "*****************************"<<endl;
        //cout << "\t In MARTINI::fragmentation mode 2: "<<endl;
        //cout << "\t\tixmax: "<<ixmax << ", izmax: "<<izmax<<endl;
        //cout << "\t\thydroZmax: "<< hydroZmax << ", hydroXmax: "<<hydroXmax <<endl;
        //cout << "\t\tNCells : "<<Ncells<< endl;
        //cout << "*****************************"<<endl;
        
        for(i=0; i < Ncells; i++)
        {
            psublist[i] = new vector<Parton>;// psublist[i] is the ith sublist that holds 
                                             // the high momentum partons that were evolved
            cvec[i] = new vector<double>;
        }
      
        countQuarkJets = 0;
        countAntiQuarkJets = 0;
      
        for ( int p_i=0; p_i<p_imax; p_i++ ) // loop over all partons in the main list 
                                             // (those above momentum threshold) 
        {                                    // to sort them into sub-lists
            id = plist->at(p_i).id();
            if (abs(id) > 22)//throw out all except q, qbar and g
                continue; 

            plist->at(p_i).status(1);
            x = plist->at(p_i).x();          // x value of position in [fm]
            y = plist->at(p_i).y();          // y value of position in [fm]
            z = plist->at(p_i).z();          // z value of position in [fm]
            t = plist->at(p_i).tFinal();     // t value of lab time in [fm/c]

            if (fabs(t) < eps)
                continue;
    
            ix = floor((hydroXmax+x)/DX+0.0001);  // x-coordinate of the cell we are in now
            iy = floor((hydroXmax+y)/DX+0.0001);  // y-coordinate of the cell we are in now
                                                  // note that x and y run from -hydroXmax to +hydroXmax
                                                  // and ix and iy from 0 to 2*hydroXmax
                                                  // hence the (hydroXmax+x or y) for both
            iz = floor((hydroZmax+z)/DZ+0.0001);// z-coordinate of the cell we are in now
  
    
            i=ix+ixmax*(iy+iz*(ixmax)); //determine sublist number by position
            if (i > Ncells-1)
            {
                cout << "*****************************"<<endl;
                cout << "x=" << x << ", y=" << y << ", z=" << z << ", t=" << t << endl;
                cout << "ix=" << ix << ", iy=" << iy << ", iz=" << iz << ", i=" << i << " out of " << Ncells << endl;
                cout << "\t In MARTINI::fragmentation mode 2: "<<endl;
                cout << "\t\tixmax: "<<ixmax << ", izmax: "<<izmax<<endl;
                cout << "\t\thydroZmax: "<< hydroZmax << ", hydroXmax: "<<hydroXmax <<endl;
                cout << "\t\tNCells : "<<Ncells<< endl;
                cout << "*****************************"<<endl;
                exit(-1);
            }
            if (plist->at(p_i).p().pAbs()<2.) // dont add parton under 2 GeV
                continue; 
            
            if (abs(id)<21)
            {
                if (id > 0) 
                    countQuarkJets++;
                else
                    countAntiQuarkJets++;
            }

            
            psublist[i]->push_back(plist->at(p_i)); // add the parton to the correct sublist.
            
            x=ix*DX-hydroXmax+DX/2.;
            y=iy*DX-hydroXmax+DX/2.;
            z=iz*DZ-hydroXmax+DZ/2.;

            if (fabs(z)>t-DZ) z=fabs(t-hydroDz)*z/fabs(z);

            cvec[i]->push_back(x); // middle of the cell
            cvec[i]->push_back(y);
            cvec[i]->push_back(z);
            cvec[i]->push_back(t);
        }

        stringCounter = 1; // gives the string a unique number

        for(i=0; i < Ncells; i++) //loop over all sublists (cells)
        {
            if (psublist[i]->size() == 0)
                continue;
            //cout << "There is/are " << psublist[i]->size() << " hard parton(s) in list " << i << endl;
            x = cvec[i]->at(0);
            y = cvec[i]->at(1);
            z = cvec[i]->at(2);
            t = cvec[i]->at(3);
            
            //cout << "x=" << x << ", y=" << y << ", z=" << z << ", t=" << t << endl;
            
            //get the temperature and flow velocity at the current position and time:
            if(fixedTemperature==0)
            {
                hydroInfo = hydroSetup->getHydroValues(x, y, z, t);
                T = hydroInfo.T;
            }
            else
                T = fixedT;
            
            if (T==0.) 
                T=0.1; // give it a minimum T
            //cout << " T in this cell =" << T << " GeV" << endl;
            
            V = DX*DX*DZ;
            
            Ngluons = (V*16./(PI*PI)*T*T*T*1.202056903/pow(hbarc,3.)); // 1.202056903 is Riemann zeta(3).
            Nquarks = (V*9./(2.*PI*PI)*T*T*T*1.202056903/pow(hbarc,3.));
            
            //cout << "N_gluons in this cell = " << Ngluons << endl;
            //cout << "N_quarks and anti-quarks in this cell = " << 6*Nquarks << endl;
            
            Igluons = round(Ngluons);
            Iquarks = round(6*Nquarks); // three flavors
            
            for (k=0; k < Ngluons; k++)
            {
                newOne.id(21);                   // add a thermal gluon
                newOne.mass(0.); 
                newOne.frozen(2);
                newOne.x(x);                     // set the new parton's initial position
                newOne.y(y);
                newOne.z(z);
                newOne.p(random->thermal(T, -1));// sample momentum from Bose distribution at temperature T
                newOne.status(0);
                psublist[i]->push_back(newOne); 
            }
            
            for (k=0; k<Nquarks; k++)
            {
                // here decide which kind of quark / antiquark it is
                rn=random->genrand64_real1();          
                if (rn<0.3333)
                {
                    id = 1;
                    mass = 0.33;
                }
                else if (rn<0.6666)
                {
                    id = 2;
                    mass = 0.33;
                }
                else 
                {
                    id = 3;
                    mass = 0.5;
                }
                rn=random->genrand64_real1();          
                if (rn<0.5) // anti-quark
                id*=-1; 
                  
                newOne.id(id); 
                newOne.mass(mass); 
                newOne.frozen(2);
                newOne.x(x);                    // set the new parton's initial position
                newOne.y(y);
                newOne.z(z);
                newOne.p(random->thermal(T, 1));// sample momentum from Fermi distribution at temperature T
                newOne.status(0);
                psublist[i]->push_back(newOne); 
            }
            //  shuffle order (so we don't always start with gluons but treat every parton on equal grounds)
            for (is=0; is<psublist[i]->size(); is++) 
            {
                rni = static_cast<int>(random->genrand64_real1()*10000) % psublist[i]->size();  // generate a random position
                temp = psublist[i]->at(is);                   // switch parton at position 'is' with that at 'rni' ...
                psublist[i]->at(is) = psublist[i]->at(rni); 
                psublist[i]->at(rni) = temp;
            }
         
            //is = 0; 
            //while (psublist[i]->begin()+is<psublist[i]->end())  // remove heavy quarks (keep photons)
            //{
            //  if (psublist[i]->at(is).id()>22) 
            //  {
            //      psublist[i]->erase(psublist[i]->begin()+is);
            //      is--;
            //  }
            //  is++;
            //}

            for (k=0; k < psublist[i]->size(); k++)  // loop over all partons in the sublist and reset color links
            {
                psublist[i]->at(k).col(0);
                psublist[i]->at(k).acol(0);
                //cout << "id=" << psublist[i]->at(k).id() << endl;
            }

            stringCounter++;
      
            for (int tries = 0 ; tries<4; tries++) // try to connect among thermal partons (only if there aren't enough add one later)
            {
                for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist
                {
                    //jets have status 1, medium partons have status 0
                    if (psublist[i]->at(k).status()==1) //  it is a jet
                    {
                        if (psublist[i]->at(k).id()==21) //   it is a gluon
                        {
                            if (psublist[i]->at(k).acol()!=0 && psublist[i]->at(k).col()!=0) continue; // is already connected
                                
                            for (l=0; l<psublist[i]->size(); l++)
                            {
                                if (l==k) continue; // do not connect to itself
                                if (psublist[i]->at(l).id()!=21 && psublist[i]->at(l).id()!=22) // if quark or anti-quark
                                {
                                    if (psublist[i]->at(l).acol()!=0 || psublist[i]->at(l).col()!=0) continue; // is already connected
                                    
                                    if (psublist[i]->at(l).id()<0 && psublist[i]->at(k).col()==0) //if anti-quark
                                    {
                                        psublist[i]->at(l).acol(stringCounter); 
                                        psublist[i]->at(k).col(stringCounter);  // connect to gluon's color
                                        stringCounter++;
                                        // cout << " 1 doing jg and a-q" << endl;
                                    }
                                    else if(psublist[i]->at(l).id()>0 && psublist[i]->at(k).acol()==0)// if quark
                                    {
                                        psublist[i]->at(l).col(stringCounter);
                                        psublist[i]->at(k).acol(stringCounter); // connect to gluon's anti-color
                                        stringCounter++;
                                        //cout << " 2 doing jg and q" << endl;
                                    }
                                }
                                else // if gluon 
                                {
                                    if (psublist[i]->at(l).acol()!=0 && psublist[i]->at(l).col()!=0) continue; // is already connected
                                    if(psublist[i]->at(k).status()>=0) // for jets and thermal (unnecessary test)
                                    {
                                        rn=random->genrand64_real1();    
                                        if (rn<0.3333) // connect color end
                                        {
                                            if (psublist[i]->at(l).acol()==0 && psublist[i]->at(k).col()==0)
                                            {
                                                psublist[i]->at(k).col(stringCounter);
                                                psublist[i]->at(l).acol(stringCounter);
                                                stringCounter++;
                                                //cout << " 3 doing jg and g color end" << endl;
                                            }
                                        }
                                        else if (rn<0.6666) // connect anti-color end
                                        {
                                            if (psublist[i]->at(l).col()==0 && psublist[i]->at(k).acol()==0)
                                            {
                                                psublist[i]->at(k).acol(stringCounter); 
                                                psublist[i]->at(l).col(stringCounter);
                                                stringCounter++;
                                                //cout << " 4 doing jg and g anti-color end" << endl;
                                            }
                                        }
                                        else //connect both ends if possible
                                        {
                                            if (psublist[i]->at(l).acol()==0 && psublist[i]->at(k).col()==0)
                                            {
                                                psublist[i]->at(k).col(stringCounter);
                                                psublist[i]->at(l).acol(stringCounter);
                                                stringCounter++;
                                                //cout << " 5 doing jg and g color end when both" << endl;
                                            }
                                            if (psublist[i]->at(l).col()==0 && psublist[i]->at(k).acol()==0)
                                            {
                                                psublist[i]->at(k).acol(stringCounter); 
                                                psublist[i]->at(l).col(stringCounter);
                                                stringCounter++;
                                                //cout << " 6 doing jg and g anti-color end when both" << endl;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else // it is a jet quark or anti-quark
                        {
                            if (psublist[i]->at(k).acol()!=0 || psublist[i]->at(k).col()!=0) continue; // is already connected
                            for (l=0; l<psublist[i]->size(); l++)
                            {
                                if (l==k) continue; // do not connect to itself
                                if (psublist[i]->at(l).id()!=21 && psublist[i]->at(l).id()!=22) // if quark or anti-quark
                                {
                                    if (psublist[i]->at(l).acol()!=0 || psublist[i]->at(l).col()!=0) continue; // is already connected
                                    if (psublist[i]->at(k).id()>0 && psublist[i]->at(l).id()<0 && psublist[i]->at(k).col()==0)
                                    //if quark and thermal anti-quark
                                    {
                                        psublist[i]->at(l).acol(stringCounter); 
                                        psublist[i]->at(k).col(stringCounter);  
                                        stringCounter++;
                                        //cout << " 7 doing jq and a-q" << endl;
                                    }
                                    else if (psublist[i]->at(k).id()<0 && psublist[i]->at(l).id()>0 && psublist[i]->at(k).acol()==0) 
                                    //if anti-quark and thermal quark
                                    {
                                        psublist[i]->at(l).col(stringCounter);
                                        psublist[i]->at(k).acol(stringCounter);
                                        stringCounter++;
                                        //cout << " 8 doing ja-q and q" << endl;
                                    }
                                }
                                else // if thermal gluon
                                {
                                    if (psublist[i]->at(k).id()>0 && psublist[i]->at(l).acol()!=0) continue;// is already connected
                                    if (psublist[i]->at(k).id()<0 && psublist[i]->at(l).col()!=0) continue; // is already connected
                                    
                                    if (psublist[i]->at(k).id()>0 && psublist[i]->at(k).col()==0) // if jet quark
                                    {
                                        psublist[i]->at(k).col(stringCounter);
                                        psublist[i]->at(l).acol(stringCounter);
                                        stringCounter++;
                                        //cout << " 9 doing jq and g" << endl;
                                    }
                                    else if (psublist[i]->at(k).id()<0 && psublist[i]->at(k).acol()==0) // if jet anti-quark
                                    {
                                        psublist[i]->at(k).acol(stringCounter);
                                        psublist[i]->at(l).col(stringCounter);
                                        //cout << " 10 doing ja-q and g" << endl;
                                        stringCounter++;
                                    }
                                }
                            }
                        }
                    }
                } // end loop over all partons in the sublist
            } // end loop over tries

            // finally connect open ends of gluons - careful here

            for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist
            {
                if (psublist[i]->at(k).id()==21 and 
                    ((psublist[i]->at(k).col()!=0 and psublist[i]->at(k).acol()==0) or 
                    (psublist[i]->at(k).acol()!=0 and psublist[i]->at(k).col()==0)) )
                {
                    if (psublist[i]->at(k).status()==1)
                      //cout << " gluon jet needs extra string attached!" << endl;
                    rn=random->genrand64_real1();         // add a thermal quark or anti-quark
                    if (rn<0.4)
                    {
                        id = 1;
                        mass = 0.33;
                    }
                    else if (rn<0.8)
                    {
                        id = 2;
                        mass = 0.33;
                    }
                    else 
                    {
                        id = 3;
                        mass = 0.5;
                    }
                    
                    newOne.mass(mass); 
                    newOne.frozen(2);
                    newOne.x(x);                              // set the new parton's initial position
                    newOne.y(y);
                    newOne.z(z);
                    newOne.p(random->thermal(T, 1)); //T         // sample momentum from Fermi distribution at temperature T
                    newOne.status(0);
                    
                    if (psublist[i]->at(k).acol()==0) // add a thermal quark and connect to it
                    {
                        //cout << "adding id " << id << " parton" << endl;
                        newOne.id(id); 
                        newOne.col(stringCounter);
                        newOne.acol(0);
                        psublist[i]->at(k).acol(stringCounter);
                        stringCounter++;
                        psublist[i]->push_back(newOne); 
                    }

                    if (psublist[i]->at(k).col()==0) // add a thermal anti-quark and connect to it
                    {
                        //cout << "adding id " << -id << " parton" << endl;
                        newOne.id(-id); 
                        newOne.acol(stringCounter);
                        newOne.col(0);
                        psublist[i]->at(k).col(stringCounter);
                        stringCounter++;
                        psublist[i]->push_back(newOne); 
                    }
                }
            }
            is = 0;
            while (psublist[i]->begin()+is<psublist[i]->end())  // remove unconnected thermal partons
            {
                if (psublist[i]->at(is).id()!=22 and psublist[i]->at(is).col()==0 and psublist[i]->at(is).acol()==0)
                {
                    //cout << "erasing position " << is << " with id=" << psublist[i]->at(is).id() << " and col/acol=" 
                    //   << psublist[i]->at(is).col() << "/" << psublist[i]->at(is).acol() << endl;
                    psublist[i]->erase(psublist[i]->begin()+is);
                    is--;
                }
                is++;
            }

            //for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist and print content
            //{
            //  if (psublist[i]->at(k).status() == 1)
            //      cout << " cell " << i << ": jet " << k << " has ID " << psublist[i]->at(k).id() 
            //           << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
            //  else if (psublist[i]->at(k).status() == 0)
            //      cout << " cell " << i << ": thermal parton " << k << " has ID " << psublist[i]->at(k).id() 
            //           << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
            //  else
            //      cout << " cell " << i << ": weird parton " << k << " has ID " << psublist[i]->at(k).id() 
            //           << " and has col=" << psublist[i]->at(k).col() << " and acol=" << psublist[i]->at(k).acol() << endl;
            //}

            for (k=0; k<psublist[i]->size(); k++)   // loop over all partons in the sublist and append to pythia list
            {
                status = psublist[i]->at(k).status();
                if (status<1 || status>100) status = 1; // if status was not set, cure this here
                //if (id==22){
                //cout << "status=" << status << endl;
                //check++;
                //}
                id = psublist[i]->at(k).id();
                if (abs(id)<4)
                {
                    if (id<0) countaq++;
                    if (id>0) countq++;
                    if (id==3) countStrange++;
                    if (id==-3) countAStrange++;
                }
                pvec = psublist[i]->at(k).p();
                col = psublist[i]->at(k).col();
                acol = psublist[i]->at(k).acol();
                mass = psublist[i]->at(k).mass();
                //cout << k << " ID=" << psublist[i]->at(k).id() 
                //   << " col=" << col << " acol=" << acol << " mass=" << mass << " p=" << pvec << endl;
                pythia.event.append(id,status,col,acol,pvec,mass); // copy existing parton into event record
            }
        } // end loop over all sublists


        //if(check>0) pythia.event.list(); 
        if(pythia.forceHadronLevel()) numEvents+=1; // does the fragmentation and the rest (if success add 1 to total events)
        //if(check>0) pythia.event.list(); 
        //cout << "q-jets " << countQuarkJets << endl;
        //cout << "qbar-jets " << countAntiQuarkJets << endl;
        //cout << "s=" << countStrange << endl;
        //cout << "sbar=" << countAStrange << endl;
        //cout << "q's=" << countq << endl;
        //cout << "qbar's=" << countaq << endl;
        
        for(i=0; i<Ncells; i++)
        {
            delete psublist[i];    // psublist[i] is the ith sublist that holds the high momentum partons that were evolved
            delete cvec[i];
        }
        delete psublist;
        delete cvec;
        
        return numEvents;
    }
}

bool MARTINI::readFile(string fileName, bool warn) 
{
  // Open file with updates.
  const char* cstring = fileName.c_str();
  ifstream is(cstring);  
  if (!is) {
    info.errorMsg("Error in MARTINI::readFile: did not find file", fileName);
    return false;
  }
  // Read in one line at a time.
  string line;
  bool accepted = true;
  while ( getline(is, line) ) 
    {
      // Process the line.
      if ( !readString( line, false ) ) // if it isn't found do not warn but pass on to PYTHIA - PYTHIA can warn if wanted
    if (!pythia.readString( line, warn ))
      accepted = false;
      // Reached end of input file.
    };
  return accepted;
}

bool MARTINI::readString(string line, bool warn) 
{
  // If empty line then done.
  if (line.find_first_not_of(" ") == string::npos) return true;

  // If first character is not a letter/digit, then taken to be a comment.
  int firstChar = line.find_first_not_of(" ");
  if (!isalnum(line[firstChar])) return true; 

  // Send on to Settings.
  return settings.readString(line, warn);
}

void MARTINI::initPythia()
{
    // In 8.3 the beams/CM energy are set via readString BEFORE init(),
    // and init() takes no arguments. The old pythia.init(2212,2212,cmEnergy),
    // setJetpTmin/max, and the *CrossSection() getters were all removed/patched.
    stringstream eCMstr;
    eCMstr << "Beams:eCM = " << cmEnergy;
    pythia.readString("Beams:idA = 2212");   // beam is a proton; the nuclear +
    pythia.readString("Beams:idB = 2212");   // isospin content comes from the EPPS21 grid
    pythia.readString(eCMstr.str());

    if ( nuclearEffects != 0 )
    {
        int Zn = glauber->get_Projetile_Z();
        int An = glauber->get_Projetile_A();
        auto pdfA = std::make_shared<NuclearAveragePDF>(PDFname, PDFmember, Zn, An);
        auto pdfB = std::make_shared<NuclearAveragePDF>(PDFname, PDFmember, Zn, An);
        pythia.setPDFPtr(pdfA, pdfB);
        cout << "[MARTINI]: nuclear PDF " << PDFname << " member " << PDFmember
             << ", average nucleon with Z=" << Zn << " A=" << An << endl;
    }

    pythia.init();   // MUST be the final PYTHIA configuration call
    cout << endl << "PYTHIA initialized with sqrt(s_NN)=" << cmEnergy << " GeV." << endl;
}


//void MARTINI::initPythia()
//{
//    // let PYTHIA know my desired infrared cutoff on the jet momentum for the calculation of the jet cross section
//    pythia.setJetpTmin(jetpTmin);
//    pythia.setJetpTmax(jetpTmax);
//      
//    // initialize PYTHIA in the CM frame.
//    pythia.init( 2212, 2212, cmEnergy ); // 2 protons, (if 14000 each 7000 GeV) /2212 is proton
//    cout << endl << "PYTHIA initialized by MARTINI with sqrt(s)=" << cmEnergy << " GeV." << endl;
//    
//    // set the relevant cross sections that PYTHIA was so kind to compute
//    elasticXSec = pythia.elasticCrossSection();
//    totalXSec = pythia.totalCrossSection();
//    jetXSec = pythia.jetCrossSection();
//    inelasticXSec = totalXSec-elasticXSec;
//    
//    // output the cross sections
//    cout << endl << "[MARTINI::initPythia]:" << endl;
//    cout << "Jet cross section = " << scientific << setprecision(6) << setw(12) << jetXSec << " mb, with a pTmin = " << jetpTmin << " GeV." << endl;
//    cout << "Total cross section = " << totalXSec << " mb." << endl;
//    cout << "Total inelastic cross section = " << inelasticXSec << " mb." << endl;
//  
//}

bool MARTINI::init(int path_number)
{
    /// set all parameters as they were initialized
    int seed = 0;
    if (path_number > 0)
    {
        seed += path_number;
    }
    // For FIC:
    file_number = seed%5 +1;
    cout << "file_number = " << file_number << endl;
    seed *= 10000;
    pCut = settings.parm("General:pCut");
    eLossCut = settings.parm("General:eLossCut");
    fixedT = settings.parm("General:Temperature");
    maxTime = settings.parm("General:MaxTime");
    dtfm = settings.parm("General:TimeStep");
    alpha_s = settings.parm("General:AlphaS");
    runs = settings.mode("General:Events");
    Nf = settings.mode("General:NumberOfFlavors");
    jetpTmin = settings.parm("General:JetPTMin");
    jetpTmax = settings.parm("General:JetPTMax");
    fullEvent = settings.mode("General:FullEvent");
    moveBeforeTau0 = settings.mode("General:MoveBeforeTau0");
    cmEnergy = settings.parm("General:cmEnergy");
    fullVacuumShower = settings.flag("General:FullVacuumShower");
    Ldependence = settings.mode("General:Ldependence");
    fragmentationMethod = settings.mode("General:FragmentationMethod");
    rateSelector = settings.mode("General:RadiativeRateSet");
    examineHQ = settings.mode("General:examineHQ");
    cout << "examineHQ = " << examineHQ << endl;
    HQ_rate_type = settings.mode("General:HQ_rate_type");
    HQ_Lattice_rate_factor = settings.parm("General:HQ_Lattice_rate_factor");
    HQ_scat_rad_flag = settings.mode("General:HQ_scat_rad_flag");
    HQ_rate_path = settings.word("General:HQ_rate_path");
    coal_prob_path = settings.word("General:coal_prob_path");
    do_coalesence = settings.flag("General:do_coalesence");
    coal_sigma = settings.parm("General:coal_sigma");
    coal_N_factor = settings.parm("General:coal_N_factor");
    cout << "HQ_rate_type = " << HQ_rate_type << endl;
    if (HQ_rate_type == 1 || HQ_rate_type ==2)
    {
      if (HQ_scat_rad_flag == 0)
      {
        cout << "Including both scattering and radiation" << endl;
      } else if (HQ_scat_rad_flag == 1)
      {
        cout << "Including only scattering" << endl;
      } else if (HQ_scat_rad_flag == 2)
      {
        cout << "Including only radiation" << endl;
      } else
      {
        cout << "Invalid option of HQ_scat_rad_flag " << HQ_scat_rad_flag;
        exit(1);
      }
    }
    setInitialPToZero = settings.mode("General:setInitialPToZero");
    cout << "setInitialPToZero = " << setInitialPToZero << endl;
    charmWidth = settings.parm("General:charmWidth");
    bottomWidth = settings.parm("General:bottomWidth");
    T_C_HQ = settings.parm("General:T_C_HQ");
    TwoPiTD_HQ = settings.parm("General:TwoPiTD_HQ");
    totalHQXSec = settings.parm("General:totalHQXSec");
    tauEtaCoordinates = settings.mode("General:tauEtaCoordinates");
    cout << "tauEtaCoordinates = " << tauEtaCoordinates << endl;

    /***** Jet Conversion Additions: Rouz *****/
    useRateTableforConversion = settings.flag("General:UseConvPhotonsTable");
    cout << "General:UseConvPhotonsTable = "<<useRateTableforConversion<<endl;
    if (useRateTableforConversion && (alpha_s > 0.32 || alpha_s < 0.28) )
    {
      cout<<"Jet Conversion using interpolation of the rate table is only available for alpha_s=0.3. Exiting." << endl;
      exit(-1);
    }
    if (useRateTableforConversion && alpha_s > 0.29 && alpha_s < 0.31)
    {
      /* I put alpha_s in (0.29, 0.31) range because we usually take alpha_s at RHIC or LHC energies to be
      * one of the following values: 0.2, 0.27, 0.3, 0.4 and I didn't want to perform an exact equality 
      * test with a double.
      * Now add the codes to read the files here. There are two files to be read
      *   1. grid information file, containing grid start and end points in 3 dimensions & grid spacing
      *   2. grid file itself.
      * Start by reading in the grid information:
      */
     cout<<"Reading conversion table data: "<<endl;
     string table_settings = settings.word("General:conversionTableSettings");
     string grid_of_rates = settings.word("General:conversionTableFile");
     cout<<" Table Settings at: "<<table_settings<<endl;
     cout<<" grid at: "<<grid_of_rates<<endl;
     angularConvGrid = new ConversionAngularGrid(table_settings, grid_of_rates);
     cout<<"Read in the jet conversion grid file."<<endl;
    }
    // *************** Done with conversion photon table stuff **********************// 
    // Check if we want a perturbative calculation:
    pertCalc = settings.flag("General:PerturbativeCalculation");
    exaggerate_brem = settings.parm("General:IncreaseRateBrem");
    exaggerate_conv = settings.parm("General:IncreaseRateConv");
    cout<<"General:PerturbativeCalculation: "<<pertCalc<<endl;
    //if ( pertCalc )
    //{
        //// perturbative calculation
        //pTmin_pert  = settings.parm("General:pTmin");
        //pTmax_pert  = settings.parm("General:pTmax"); 
        //nbins_pert  = settings.parm("General:NumBins"); 
        //etaCut_pert = settings.parm("General:etaCut");
        //conversion_photons = new TH1D("conv_gamma_hist", "pT hist of conv gamma", nbins_pert, pTmin_pert, pTmax_pert);
        //amy_photons        = new TH1D("amy_gamma_hist", "pT hist of amy gamma", nbins_pert, pTmin_pert, pTmax_pert);

        //TH1::SetDefaultSumw2(); // so that ROOT would keep stat errors
    //}
    //-----------------------------------------------------------------------------//
    if (rateSelector < 0 || rateSelector > 4)
    {
        cout << "The chosen radiative rates are not available at the moment. Exiting." << endl;
        exit(1);
    }
    if (settings.flag("General:RadiativeProcesses")) doRadiative=1;
    else doRadiative=0;
    if (settings.flag("General:ElasticCollisions")) doElastic=1;
    else doElastic=0;
    if (settings.flag("General:TransferTransverseMomentum")) transferTransverseMomentum=1;
    else transferTransverseMomentum=0;
    if (settings.flag("General:Evolution")) evolution=1;
    else evolution=0;
    if (settings.flag("General:EvolutionShear")) evolutionshear=1;
    else evolutionshear=0;
    if (settings.flag("General:EvolutionBulk")) evolutionbulk=1;
    else evolutionbulk=0;
    if (settings.flag("General:FixedEnergy")) fixedEnergy=1;
    else fixedEnergy=0;
    if (settings.flag("General:FixedTemperature")) fixedTemperature=1;
    else fixedTemperature=0;
    if (settings.flag("General:Fragmentation")) fragmentationSwitch=1;
    else fragmentationSwitch=0;
    if (settings.flag("General:PhotonProduction")) photonSwitch=1;
    else photonSwitch=0;
    if (settings.flag("General:prehydroevolution")) prehydroevolution=1;
    else prehydroevolution=0;
    if (settings.flag("PDF:NuclearEffects")) nuclearEffects=1;
    else nuclearEffects=0;

    useLHAPDF = settings.flag("PDF:useLHAPDF");
    getRecoil = settings.flag("General:GetRecoil"); 
    cout << "getRecoil = " << getRecoil << endl;
    recoilCut = settings.parm("General:RecoilCut");
    cout << "recoilCut = " << recoilCut << endl;
    outputSource = settings.flag("General:OutputSource"); 
    cout << "outputSource = " << outputSource << endl;
    trackHistory = settings.flag("General:TrackHistory"); 
    initialXjet = settings.parm("General:InitialXjet");
    initialYjet = settings.parm("General:InitialYjet");
    initialPXjet = settings.parm("General:InitialPXjet");
    initialPYjet = settings.parm("General:InitialPYjet");
    allFromCenter = settings.flag("General:allFromCenter");
    pCutPropToT = settings.flag("General:pCutPropToT");
    formationTime = settings.flag("General:FormationTime");
    cout << "formationTime = " << formationTime << endl;
    runningElas = settings.flag("General:RunningElas");
    cout << "runningElas = " << runningElas << endl;
    runningRad = settings.flag("General:RunningRad");
    cout << "runningRad = " << runningRad << endl;
    renormFacRad = settings.parm("General:RenormFacRad");
    cout << "renormFacRad = " << renormFacRad << endl;
    renormFacElas = settings.parm("General:RenormFacElas");
    cout << "renormFacElas = " << renormFacElas << endl;
    prehydroTau0 = settings.parm("PreHydro:Tau0");
    prehydroDtau = settings.parm("PreHydro:Dtau");
    prehydroXmax = settings.parm("PreHydro:Xmax");
    prehydroDx = settings.parm("PreHydro:Dx");
    hydroTau0 = settings.parm("Hydro:tau0");
    hydroTauMax = settings.parm("Hydro:taumax");
    hydroDtau = settings.parm("Hydro:dtau");
    hydroDx = settings.parm("Hydro:dx");
    //Hydrodynamical output saved in tau-eta coordinates will simply use hydroDz and hydroZmax:
    hydroDz = settings.parm("Hydro:dz");
    hydroXmax = settings.parm("Hydro:xmax");
    hydroZmax = settings.parm("Hydro:Zmax");
    hydro_nskip_tau = settings.parm("Hydro:nskip_tau");
    hydro_nskip_x = settings.parm("Hydro:nskip_x");   // for furture
    hydro_nskip_z = settings.parm("Hydro:nskip_z");   // for furture
    hydroTfinal = settings.parm("Hydro:Tfinal");
    hydroWhichHydro = settings.mode("Hydro:WhichHydro");  
    hydroSubset = settings.mode("Hydro:Subset");
    hydroViscous = settings.flag("Hydro:Viscous");  
    fixedDistribution = settings.flag("Hydro:fixedDistribution");
    
    hydroDtau = hydroDtau*hydro_nskip_tau;
    hydroDx = hydroDx*hydro_nskip_x;
    hydroDz = hydroDz*hydro_nskip_z;
    
    glauberTarget = settings.word("Glauber:Target");
    glauberProjectile = settings.word("Glauber:Projectile");
    glauberImpactParam = settings.parm("Glauber:b");
    glauberIMax = settings.mode("Glauber:MaxInterpolationPoints");
    glauberEnvelope = settings.flag("Glauber:Envelope");
    
    PDFname = settings.word("PDF:LHAPDFset");
    PDFmember = settings.mode("PDF:LHAPDFmember");
    nuclearPDFname = settings.word("PDF:nuclearPDF");
    
    nbinFromFile = settings.mode("General:Nbin_from_File");
    NColl = settings.mode("General:NColl");
    binary_filename = settings.word("General:Nbin_File_Name");
    //Either the name of the evolution file, or a beginning tag for the 
    //file for the case nbinFromFile == 1. -CFY
    evolution_name = settings.word("General:evolution_name");
    pre_evolution_name = settings.word("General:pre_evolution_name");
    IPG_filename = settings.word("General:IPG_filename");
    shear_name = settings.word("General:shear_name");
    bulk_name = settings.word("General:bulk_name");
    //The number of collisions in the file, counted up from zero
    cout << "evolution_name = " << evolution_name << endl;
    cout << "binary_filename = " << binary_filename << endl;
    background_file_name = settings.word("General:Background_file_name");
    cout << "background_file_name = " << background_file_name << endl;
    Nbin = 0;
  
    if(nbinFromFile == 2 or nbinFromFile == 3 or nbinFromFile == 4)
    {
        binary_info_ptr->init(binary_filename, nbinFromFile);
        //binary_info_ptr->print_info();
        //binary_info_ptr->check_samples();
    }
  
//    if ( fixedTemperature == 0 and fixedEnergy == 0)
//        initPythia();

//    if ( fixedTemperature == 0 and fixedEnergy == 0)
        // init PYTHIA with given center of mass energy. will be changed for full AA collision.
//        initPythia();
    ///read in data files with transition rates
    if (examineHQ == 0)
    {
      import->init(rateSelector);
      rates   = new Rates(rateSelector);
    }
    if ( fixedTemperature == 0 and fixedEnergy == 0)
    {
          // initialize nuclei information
          inelasticXSec = settings.parm("General:inelasticNNXSec");
          glauber->init(inelasticXSec,glauberTarget,glauberProjectile,glauberImpactParam,glauberIMax,glauberEnvelope);
    }  

    //importLRates->init();
    
    //cout << "rate=" << importLRates->getRate(10., 5., 2.) << endl;
    
    // output of the rates into a file:
    //import->show_dGamma();
    
    /// initialize random number generator with current time as seed
    long long rnum;
    
    rnum=time(0)+seed;
    
    random->init_genrand64(rnum);
    
    // set random seed as current time
    srand(rnum);
    int now=rand()/100;
    stringstream nowstr;
    nowstr << "Random:seed = " << now;
    string timestring = nowstr.str();
    // put most usual PYTHIA setups in this if-statement
    // we only want this to be done in the production run
    // and not during a QGP brick + parton gun mode. 
    // RY June 9 2021
    if ( fixedTemperature == 0 and fixedEnergy == 0)
    {
        // PYTHIA settings:
        pythia.readString("Random:setSeed = on");

        pythia.readString(timestring);    
        //pythia.readString("Random:seed = 100");

        // finish after doing hard process
        pythia.readString("PartonLevel:all = on"); // off to only get hard process partons
        pythia.readString("HadronLevel:all = off");

        stringstream ptminstr;
        ptminstr << "PhaseSpace:pTHatMin =  " << jetpTmin;
        string ptminstring = ptminstr.str();
        pythia.readString(ptminstring);
        stringstream ptmaxstr;
        ptmaxstr << "PhaseSpace:pTHatMax =  " << jetpTmax;
        string ptmaxstring = ptmaxstr.str();
        pythia.readString(ptmaxstring);

        // very important to have this option!!
        pythia.readString("Check:event = off");
    }
    if(examineHQ == 1)
    {
        pythia.readString("HardQCD:gg2ccbar = on");
        pythia.readString("HardQCD:qqbar2ccbar = on");
        pythia.readString("HardQCD:gg2bbbar = on");
        pythia.readString("HardQCD:qqbar2bbbar = on");
    }

    if ( nuclearEffects != 0 )
    {
 //       stringstream PDFnameStr;
//        PDFnameStr << "PDF:pSet = LHAPDF6:" << PDFname;   // e.g. EPPS21nlo_CT18Anlo_208
//        pythia.readString(PDFnameStr.str());
    }
    else
    {
        if (useLHAPDF)
        {
            stringstream PDFnameStr;
            PDFnameStr << "PDF:pSet = LHAPDF6:" << PDFname;
            pythia.readString(PDFnameStr.str());
        }
        // else: PYTHIA's built-in proton PDF (NNPDF2.3 LO) is used
    }



//    // set the PDF and possible nuclear effects. note: if nuclear effects are chosen, the use of LHAPDF is enforced!
//    if ( nuclearEffects!=0 )
//    {
//        if (!useLHAPDF)
//            cout << "[MARTINI]:WARNING: using nuclear effects - turned on LHAPDF against initial settings." << endl;
//        pythia.readString("PDF:nuclearEffects = 1");
//        stringstream PDFnameStr;
//        PDFnameStr << "PDF:pSet=LHAPDF5:" << PDFname << "/" << PDFmember;
//        pythia.readString(PDFnameStr.str());
//        stringstream nuclearPDFnameStr;
//        nuclearPDFnameStr << "PDF:nuclearPDF=" << nuclearPDFname;
//        pythia.readString(nuclearPDFnameStr.str());
//        stringstream projAstr, targAstr;
//        projAstr << "PDF:proj_atomicNumber = " << glauber->get_Projetile_A();
//        targAstr << "PDF:targ_atomicNumber = " << glauber->get_Target_A();
//        pythia.readString(projAstr.str());
//        pythia.readString(targAstr.str());
//    }
//    else
//    {
//        if(useLHAPDF) 
//        {
//            stringstream PDFnameStr;
//            PDFnameStr << "PDF:pSet=LHAPDF5:" << PDFname << "/" << PDFmember;
//            pythia.readString(PDFnameStr.str());
//        }
//        pythia.readString("PDF:nuclearEffects = 0");
//    }
  
    if ( fixedTemperature == 0 and fixedEnergy == 0)
        // init PYTHIA with given center of mass energy. will be changed for full AA collision.
        initPythia();

    // read the hydro data from the data file(s)
    if (prehydroevolution == 1 and fixedTemperature==0) 
    {
        hydroSetup->readPreHydroData(prehydroTau0, hydroTau0, prehydroDtau, 
                                  prehydroXmax, prehydroDx, pre_evolution_name); 
        cout << "OK after readPreHydroData" << endl;
    }
    if (evolution and fixedTemperature==0) 
    {
        hydroSetup->readHydroData(hydroTau0, hydroTauMax, hydroDtau, 
                                  hydroXmax, hydroZmax, hydroDx, hydroDz, 
                                  hydro_nskip_tau, hydro_nskip_x, hydro_nskip_z,
                                  hydroWhichHydro, hydroTfinal, tauEtaCoordinates,
                                  evolutionshear, evolutionbulk, 
                                  evolution_name, shear_name, bulk_name); 
        cout << "OK after readHydroData" << endl;
        // update hydroTauMax according to the actual size of the hydro medium
        hydroTauMax = hydroSetup->get_hydro_tau_max();
    }
    if (examineHQ != 0)
    {
        if (HQ_rate_type == 1 || HQ_rate_type == 2)
        {
          hqRates->readHQRates(HQ_scat_rad_flag, HQ_rate_path); 
        }
	if (do_coalesence) hqRates->readCoalesenceProb(coal_prob_path);
    }
    if (trackHistory)
    {
        Tt = new vector<double>; // stores the history of the trajectory: t, x(t), y(t), z(t), dE/dt(t), dpx/dt(t), dpy/dt(t), dpz/dt(t) *for one particle*
        Tx = new vector<double>;
        Ty = new vector<double>;
        Tz = new vector<double>;
        TE = new vector<double>;
        Tpx = new vector<double>;
        Tpy = new vector<double>;
        TdEdt = new vector<double>;
        Tdpxdt = new vector<double>;
        Tdpydt = new vector<double>;
        Tdpzdt = new vector<double>;
        hydroSetup->output_temperature_evolution("temp_evo");
    }

    return true;
}

//Sangyong's addition
// things to do when there is a daughter
int MARTINI::hasDaughter(vector<Parton> *plist, int i, double T)
{
    int daughter;
    double x, y, z, rperp, kx, ky, kz, kox, koy, koz, ko;
    double kcrx, kcry, kcrz, kperp;
    
    daughter = plist->at(i).daughter();
    
    kox = plist->at(i).p_at_split().px();
    koy = plist->at(i).p_at_split().py();
    koz = plist->at(i).p_at_split().pz();
    ko = sqrt(kox*kox + koy*koy + koz*koz);
    if (ko == 0.0)
    {
        kox = plist->at(daughter).p_at_split().px();
        koy = plist->at(daughter).p_at_split().py();
        koz = plist->at(daughter).p_at_split().pz();
        ko = sqrt(kox*kox + koy*koy + koz*koz);
        if (ko == 0.0)
        {
             cout << "hasDaughter: this can'happen.\n" << endl;
             exit(0);
        }
    }
 
     x = plist->at(i).x() - plist->at(daughter).x();
     y = plist->at(i).y() - plist->at(daughter).y();
     z = plist->at(i).z() - plist->at(daughter).z();

     // k_orig cross (delta x) = k_orig cross (delta xperp)
     // | k_orig cross (delta x)| = k_orig rperp

     kcrx = y*koz - z*koy;
     kcry = z*kox - x*koz;
     kcrz = x*koy - y*kox;
     rperp = sqrt(kcrx*kcrx + kcry*kcry + kcrz*kcrz)/ko;
    
     kx = plist->at(i).p().px() - plist->at(daughter).p().px();
     ky = plist->at(i).p().py() - plist->at(daughter).p().py();
     kz = plist->at(i).p().pz() - plist->at(daughter).p().pz();
    
     // k_orig cross (delta k) = k_orig cross (delta kperp)
     // | k_orig cross (delta k)| = k_orig kperp

     kcrx = ky*koz - kz*koy;
     kcry = kz*kox - kx*koz;
     kcrz = kx*koy - ky*kox;
     kperp = sqrt(kcrx*kcrx + kcry*kcry + kcrz*kcrz)/ko;
    
     double k = plist->at(daughter).p().pAbs();
     double Crperp = 0.25*pow(k/T, 0.11);
    // two lows : 0.05 and 0.01
    // two highs: 0.08 and 2000000(if I rememeber correctly)
     if ( rperp*kperp < Crperp*hbarc ) 
     {
        return 1; // do only elastic
     }
     else  // decouple
     {
        plist->at(i).mother_of(-1); 
        plist->at(i).daughter_of(-1); 
        plist->at(daughter).mother_of(-1); 
        plist->at(daughter).daughter_of(-1); 
    
        return 0; // do rad plus elastic
     }
}// hasDaughter;


//Sangyong's addition
// things to do when there is a mother
int MARTINI::hasMother(vector<Parton> *plist, int i, double T)
{
    int mother;
    double x, y, z, rperp, kx, ky, kz, kox, koy, koz, ko;
    double kcrx, kcry, kcrz, kperp;
    
    mother= plist->at(i).mother();
    
    kox = plist->at(mother).p_at_split().px();
    koy = plist->at(mother).p_at_split().py();
    koz = plist->at(mother).p_at_split().pz();
    ko = sqrt(kox*kox + koy*koy + koz*koz);
    if( ko == 0.0 )
    {
        kox = plist->at(i).p_at_split().px();
        koy = plist->at(i).p_at_split().py();
        koz = plist->at(i).p_at_split().pz();
        ko = sqrt(kox*kox + koy*koy + koz*koz);
        if( ko == 0.0 )
        {
            cout << "hasMother: this can'happen.\n" << endl;
            exit(0);
        }
    }// if ko == 0

    x = plist->at(i).x() - plist->at(mother).x();
    y = plist->at(i).y() - plist->at(mother).y();
    z = plist->at(i).z() - plist->at(mother).z();
    
    // k_orig cross (delta x) = k_orig cross (delta xperp)
    // | k_orig cross (delta x)| = k_orig rperp
    
    kcrx = y*koz - z*koy;
    kcry = z*kox - x*koz;
    kcrz = x*koy - y*kox;
    rperp = sqrt(kcrx*kcrx + kcry*kcry + kcrz*kcrz)/ko;
    
    kx = plist->at(i).p().px() - plist->at(mother).p().px();
    ky = plist->at(i).p().py() - plist->at(mother).p().py();
    kz = plist->at(i).p().pz() - plist->at(mother).p().pz();
    
    // k_orig cross (delta k) = k_orig cross (delta kperp)
    // | k_orig cross (delta k)| = k_orig kperp
    
    kcrx = ky*koz - kz*koy;
    kcry = kz*kox - kx*koz;
    kcrz = kx*koy - ky*kox;
    kperp = sqrt(kcrx*kcrx + kcry*kcry + kcrz*kcrz)/ko;
    
    double k = plist->at(i).p().pAbs();
    double Crperp = 0.25*pow(k/T, 0.11);
    if ( rperp*kperp < Crperp*hbarc ) 
    {
        return 1; // do only elastic
    }
    else  // decouple
    {
        plist->at(i).mother_of(-1); 
        plist->at(i).daughter_of(-1); 
        plist->at(mother).mother_of(-1); 
        plist->at(mother).daughter_of(-1); 
        
        return 0; // do rad plus elastic
    }

}// hasMother;

void MARTINI::output_tracked_history(string filename)
{
    cout << "MARTINI:: output tracked histroy into file: " << filename << endl;
    ofstream history(filename.c_str());
    for(int i = 0; i < Tt->size(); i++)
    {
        history << scientific << setw(16) << setprecision(6)
                << (*Tt)[i] << "  " << (*Tx)[i] << "  " << (*Ty)[i] << "  " << (*Tz)[i] << "  "
                << (*TE)[i] << "  " << (*Tpx)[i] << "  " << (*Tpy)[i] << endl;
    }
    history.close();
}
int MARTINI::generate_pGun(vector<Parton> *plist, int gunID)
{
    // RY June 9 2021
    //cout<<"MARTINI: generate a parton gun of type "<< gunID<<endl; 
    //cout<<"Ensure fixed temperature and energy:"<<endl;
    if ( fixedTemperature == 0 or fixedEnergy == 0)
    {
        cout<<"[ERROR]: parton gun mode must be at fixed E and T. Exit."<<endl;
        exit(-1);
    }
    Parton parton;
    double mass;
    if ( gunID >= 1 and gunID <= 3)//quark
    {
        parton.col(1000000);
        parton.acol(0);
    }
    else if ( gunID <= -1 and gunID >= -3 )//anti-quark
    {
        parton.col(0);
        parton.acol(1000000);
    }
    else if ( gunID == 21 ) //gluon
    {
        parton.col(1000000);
        parton.acol(1000001);
    }
    else
    {
        string test = abs(gunID) < 1 ? "True" : "False";
        cout<<"gunID: "<<gunID<<" abs(gunID) < 1 ? "<< test<<endl;
        //it's not a u, ubar, d, dbar, s, sbar or g so
        // no: say it's not allowed and terminate.
        cout<<"[ERROR]: parton ID is not recognized. Exit."<<endl;
        exit(-1);
    }
    if ( abs(gunID) == 1 or abs(gunID) == 2) // d or u quark
        mass = 0.33;
    else if ( abs(gunID) == 3) // s quark
        mass = 0.55;
    else // gluon
        mass = 0.0;

    Vec4 pvec;
    double energy = sqrt(initialPXjet*initialPXjet + initialPYjet*initialPYjet + mass*mass);
    pvec.px(initialPXjet);
    pvec.py(initialPYjet);
    pvec.pz(0.0);
    pvec.e(energy);
    parton.mass(mass);
    parton.p(pvec);
    parton.id(gunID); 
    parton.source(10);
    parton.status(16);
    parton.x(initialXjet);
    parton.y(initialYjet);
    parton.z(0.0);
    parton.tFinal(0.0);
    parton.mother_of(-1); // mother of nobody
    parton.daughter_of(-1); // daughter of nobody
    parton.p_at_split(parton.p()); // original momentum, don't need it but might as well populate it
    parton.init_counts();
    plist->push_back(parton); 
    return 1;
}

//Begin rewrote routines for HQ by CFY -- MS

//Here, we add the routines needed for heavy quark and quarkonium production:

double MARTINI::CornellPotential(double r) {

  //In GeV:                                                                                                                                                                    
  if(r < 0.2) r = 0.2 ;
  double potential = -0.0985/r+0.812*r;
  if(potential > 1.45){
    potential = 1.45 ; }

  return potential ;
}

double MARTINI::CornellPotential(double r, int id1, int id2) {

  //In GeV:                                                                                                                                                                    
  if(r < 0.2) r = 0.2 ;
  double potential = -0.0985/r+0.812*r;

  double VCMax;
  if(fabs(id1) == 4 && fabs(id2) == 4){
    VCMax = E_CCB_BOUND-2.*CHARM_MASS;
  }
  else if(fabs(id1) == 5 && fabs(id2) == 5){
    VCMax = E_BBB_BOUND-2.*BOTTOM_MASS;
  }
  else{ VCMax = E_BBC_BOUND-CHARM_MASS-BOTTOM_MASS;}

  if(potential > VCMax){
    potential = VCMax ; }

  return potential ;
}

double MARTINI::FCornell(double r) {

  double F;
  if(r < 0.2 || r > 1.85124){ F = 0.;}
  else{ F = 0.0985/(r*r)+0.812;}

  return F;
}

double MARTINI::FCornell(double r, int id1, int id2) {

  double F;

  double VCMax;
  if(fabs(id1) == 4 && fabs(id2) == 4){
    VCMax = E_CCB_BOUND-2.*CHARM_MASS;
  }
  else if(fabs(id1) == 5 && fabs(id2) == 5){
    VCMax = E_BBB_BOUND-2.*BOTTOM_MASS;
  }
  else{ VCMax = E_BBC_BOUND-CHARM_MASS-BOTTOM_MASS;}

  double rMax = (VCMax+sqrt(VCMax*VCMax+4.*0.0985*0.812))/(2.*0.812) ;

  if(r < 0.2 || r > rMax){ F = 0.;}
  else{ F = 0.0985/(r*r)+0.812;}

  return F;
}

double MARTINI::F(double r, double T){
  double F;
  if(r < 0.2){
    F = 0.;
  }
  else{
    double c_term = 0.14657/(T/T_C_HQ-0.98) ;
    double mu_term = 0.0301+0.0063*T/T_C_HQ ;
    double exp_term = exp(-mu_term*r*r/(hbarc*hbarc)) ;
    
    //The magnitude of the force, from a fit to quenched
    //QCD data:
    F = 31.841*mu_term*exp_term
      +0.6199*exp_term/(r*r)
      +1.1258*exp_term-57.823*r*r*mu_term*exp_term
      +7.528*mu_term*r*exp_term/(T/T_C_HQ-0.98) ;
    //warning, Warning, WARNING!!!!
    F = hbarc*F;
  }
  return F;
}

void MARTINI::stepMomentumInCMFrame(double *ptot, double *pmu, double dt, double *deltax, double *pnew, double T){
  double invMass = sqrt(pmu[0]*pmu[0] -pmu[1]*pmu[1] -pmu[2]*pmu[2] -pmu[3]*pmu[3]);

  double uboost[4];
  uboost[0] = 1.;
  uboost[1] = ptot[1]/ptot[0];
  uboost[2] = ptot[2]/ptot[0];
  uboost[3] = ptot[3]/ptot[0];
  double v2 = uboost[1]*uboost[1]+uboost[2]*uboost[2]+uboost[3]*uboost[3];
  double gammaboost = 1./sqrt(1.-v2);
  for(int iu=0; iu<4; iu++){
    uboost[iu] *= gammaboost;
  }
  //The boost back:
  double uboostback[4];
  uboostback[0] = uboost[0];
  uboostback[1] = -uboost[1];
  uboostback[2] = -uboost[2];
  uboostback[3] = -uboost[3];
  
  //First, the parton's momentum in the CM frame:
  double pCMFrame[4];
  LorentzBooster(uboost, pmu, pCMFrame);
  
  //Next, the timestep estimated in this frame:
  double CMdt = dt*sqrt(ptot[0]*ptot[0]-ptot[1]*ptot[1]-ptot[2]*ptot[2]-ptot[3]*ptot[3] )/ptot[0];
  //Now, in the CM frame:
  double CMdeltax[4];
  LorentzBooster(uboost, deltax, CMdeltax);

  //The norm of \vec CMdeltax is also needed:
  double r = sqrt(CMdeltax[1]*CMdeltax[1]+CMdeltax[2]*CMdeltax[2]+CMdeltax[3]*CMdeltax[3] );

  double force;
  if(T > T_C_HQ){
    force = F(r, T);
  }
  else{force = FCornell(r);}

  double F_x, F_y, F_z;
  if(force == 0.){F_x = F_y = F_z = 0.;}
  else{
    F_x = force*CMdeltax[1]/r;
    F_y = force*CMdeltax[2]/r;
    F_z = force*CMdeltax[3]/r;
  }

  pCMFrame[1] += F_x*CMdt*hbarc;
  pCMFrame[2] += F_y*CMdt*hbarc;
  pCMFrame[3] += F_z*CMdt*hbarc;
  pCMFrame[0] = sqrt(invMass*invMass +pCMFrame[1]*pCMFrame[1] +pCMFrame[2]*pCMFrame[2] +pCMFrame[3]*pCMFrame[3]);
  
  LorentzBooster(uboostback, pCMFrame, pnew);

  return;
}

double MARTINI::rInCMFrame(double *p1, double *p2, double *deltaX){
  double ptot[4];
  for(int ip=0; ip<4; ip++){
    ptot[ip] = p1[ip]+p2[ip];
  }

  double uboost[4];
  uboost[0] = 1.;
  uboost[1] = ptot[1]/ptot[0];
  uboost[2] = ptot[2]/ptot[0];
  uboost[3] = ptot[3]/ptot[0];
  double gammaboost = 1./sqrt(1.-uboost[1]*uboost[1]-uboost[2]*uboost[2]-uboost[3]*uboost[3]);
  for(int iu=0; iu<4; iu++){
    uboost[iu] *= gammaboost;
  }

  double CMDeltaX[4];
  LorentzBooster(uboost, deltaX, CMDeltaX);

  double r = sqrt( CMDeltaX[1]*CMDeltaX[1] + CMDeltaX[2]*CMDeltaX[2] + CMDeltaX[3]*CMDeltaX[3] );

  return r;
}

int MARTINI::generateEventHeavyQuarks(vector<Parton> *plist)
{
  pythia.event.reset();
  if (fullEvent==1)
    {
      // Sample Nu (which gives the number of nucleons at radial position r) for both nuclei
      // then lay them on top of each other and sample the number of charm and bottom events with 
      // area = totalHQXSec.
      
      //cout << "Number of binary collisions = " << glauber->TAB() << endl; 
      //cout << " A=" << glauber->nucleusA() << endl;

      Parton parton;                               // new Parton object
      double b = glauberImpactParam;
      int posx;
      int posy;
      int n1=0;
      int n2=0;
      double r = 1.2*pow(glauber->get_Projetile_A(), 1./3.);
//      double r = 1.2*pow(glauber->nucleusA(),1./3.);

      nucleusA.clear();
      nucleusB.clear();
      sampleTA();                                  // populate the lists nucleusA and nucleusB with position data of the nucleons

      double cellLength = sqrt(inelasticXSec*0.1); // compute cell length in fm (1 mb = 0.1 fm^2)
      double latticeLength = 4.*r;                 // spread the lattice 2r in both (+/-) directions (total length = 4r)
      int ixmax = ceil(latticeLength/cellLength);
      ixmax*=2;
      double xmin = -latticeLength;
      
      int nucALat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus A
      int nucBLat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus B
      int collLat[ixmax][ixmax];                   // lattice with the positions of the interactions -> for information only

      int countCollisionsA[100];
      int countCollisionsB[100];

      for (int i=0; i<100; i++)
  {
    countCollisionsA[i]=0;
    countCollisionsB[i]=0;
  }

      for (int i = 0; i < ixmax; i++)              // initialize cells: zero nucleons in every cell
  for (int j = 0; j < ixmax; j++)
    {
      nucALat[i][j] = 0;
      nucBLat[i][j] = 0;
      collLat[i][j] = 0;
    }

      for (int i = 0; i < nucleusA.size(); i++)    // fill nucleons into cells
  {
    posx = floor((nucleusA.at(i).x-xmin-b/2.)/cellLength);
    posy = floor((nucleusA.at(i).y-xmin)/cellLength);
    nucALat[posx][posy] += 1;
        
    posx = floor((nucleusB.at(i).x-xmin+b/2.)/cellLength);
    posy = floor((nucleusB.at(i).y-xmin)/cellLength);
    nucBLat[posx][posy] += 1;

  }

      int done;
      int numberOfpp;
      int eventNumber=0;
      numberOfpp=0;
      //cout << "jetXSec=" << jetXSec << endl;
      //cout << "inelasticXSec=" << inelasticXSec << endl;
      int Projectile_A=glauber->get_Projetile_A();
      int Target_A=glauber->get_Projetile_A();
      int Projectile_Z = glauber->get_Projetile_Z();
      int Target_Z = glauber->get_Projetile_Z();
      for (int ix = 0; ix < ixmax; ix++) 
  for (int iy = 0; iy < ixmax; iy++)
    {
      for (int i = 0; i < nucALat[ix][iy]; i++)
        for (int j = 0; j < nucBLat[ix][iy]; j++)
    {
      //cout << "ix+b/2=" << static_cast<int>(ix+b/2.) << " # there=" << nucALat[static_cast<int>(ix+b/2.)][iy] << endl;
      if ( random->genrand64_real1() < totalHQXSec/inelasticXSec )
        {
          countCollisionsA[i]+=1;
          countCollisionsB[j]+=1;
          double xPositionInCell = (random->genrand64_real1())*cellLength;
          double yPositionInCell = (random->genrand64_real1())*cellLength;
          // see if colliding nucleons are neutrons
          if ( random->genrand64_real1() > (double)Projectile_Z/(double)Projectile_A ) n1=2112;
          else n1=2212;
          if ( random->genrand64_real1() > (double)Target_Z/(double)Target_A ) n2=2112;
          else n2=2212;

          int pythiaWorked = 0;
          while(pythiaWorked == 0){
          pythia.next();
//      pythiaWorked = pythia.next(n1,n2);
          }                                                // generate event with pythia

          done = 0;
          for (int ip = 0; ip < pythia.event.size(); ++ip) 
      {
        if (pythia.event[ip].status()>0 && (abs(pythia.event[ip].id())==4 || abs(pythia.event[ip].id()) == 5) )
          // If the parton is final, and is a heavy quark, then put it into the list.
          {
            parton.id(pythia.event[ip].id());                     // set parton id
            parton.status(pythia.event[ip].status());             // set parton status
            parton.mass(pythia.event[ip].m());                // set mass
            if (fixedTemperature==0)
        {
          //Random numbers added to the initial positions:
          double dx, dy, dz;
          dx = dy = dz = 0.;
          if(abs(parton.id()) == 4 && charmWidth > 0.01){
            dx = gsl_ran_gaussian(gsl_rand, charmWidth);
            dy = gsl_ran_gaussian(gsl_rand, charmWidth);
          }
          if(abs(parton.id()) == 5 && bottomWidth > 0.01){
            dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
            dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
          }

          if( moveBeforeTau0 == 1){
            Vec4 partonP = pythia.event[ip].p();
            double pzP = partonP.pz();
            double EP = sqrt(pythia.event[ip].m()*pythia.event[ip].m()
                 +partonP.px()*partonP.px()
                 +partonP.py()*partonP.py()
                 +partonP.pz()*partonP.pz() );
            double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
            dx += partonP.px()*dt/EP;
            dy += partonP.py()*dt/EP;
            dz += partonP.pz()*dt/EP;
          }

          parton.x(xmin+(ix)*cellLength+xPositionInCell+dx);       // set position
          parton.y(xmin+(iy)*cellLength+yPositionInCell+dy);
//          parton.xini(xmin+(ix)*cellLength+xPositionInCell+dx);       // set position
//          parton.yini(xmin+(iy)*cellLength+yPositionInCell+dy);
          parton.z(dz);
//          parton.zini(dz);
//          parton.tini(hydroTau0);
          //Set the old positions to be the current ones, for the first timestep:
          parton.xOld(xmin+(ix)*cellLength+xPositionInCell+dx);       // set position
          parton.yOld(xmin+(iy)*cellLength+yPositionInCell+dy);
          parton.zOld(dz);
          if ( done == 0 )
            {
              posx = floor((parton.x()-xmin)/cellLength);
              posy = floor((parton.y()-xmin)/cellLength);
              collLat[posx][posy] += 1;
                          
              //cout << numberOfpp << endl;
              numberOfpp++;
              done = 1;
            }
        }
            else
        {
          double dx, dy, dz;
          dx = dy = dz = 0.;
          if(abs(parton.id()) == 4){
            dx = gsl_ran_gaussian(gsl_rand, charmWidth);
            dy = gsl_ran_gaussian(gsl_rand, charmWidth);
          }
          if(abs(parton.id()) == 5){
            dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
            dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
          }

          if( moveBeforeTau0 == 1){
            Vec4 partonP = pythia.event[ip].p();
            double pzP = partonP.pz();
            double EP = sqrt(pythia.event[ip].m()*pythia.event[ip].m()
                 +partonP.px()*partonP.px()
                 +partonP.py()*partonP.py()
                 +partonP.pz()*partonP.pz() );
            double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
            dx += partonP.px()*dt/EP;
            dy += partonP.py()*dt/EP;
            dz += partonP.pz()*dt/EP;
          }

          parton.x(dx);     
          parton.y(dy);
//          parton.xini(dx);     
//          parton.yini(dy);
          parton.z(dz);
//          parton.zini(dz);
//          parton.tini(0.);
          //Set the old positions to be the current ones, for the first timestep:
          parton.xOld(dx);     
          parton.yOld(dy);
          parton.zOld(dz);
        }
            parton.col(pythia.event[ip].col());                 // set color 
            parton.acol(pythia.event[ip].acol());                 // set anti-color
            parton.frozen(0);                                     // parton is not frozen (will evolve)
            parton.formed(0);                                     // parton is not formed immediately (will evolve)
            parton.inhydro(0);                                     // parton is not formed immediately (will evolve)
            if(setInitialPToZero == 0){
        parton.p(pythia.event[ip].p());  // set momentum
        parton.Initialp(pythia.event[ip].p());  // set momentum
        parton.pOld(pythia.event[ip].p());  // set momentum
//        parton.pini(pythia.event[ip].p());                    // set initial momentum
            }
            else{
        parton.p(0.,0.,0.);
        parton.Initialp(0.,0.,0.);
        parton.pOld(0.,0.,0.);
//        Vec4 pinitial(pythia.event[ip].mass(),0.,0.,0.);
//        parton.pini(pinitial);                    // set initial momentum
            }
            parton.tFinal(0.);                                    // Set the initial freezeout time to zero
//            parton.eventNumber(eventNumber);
//            parton.elasticCollisions(0);                          // set initial no. of el colls to zero
            parton.source(0);                                     // All initial partons have itsSource=0
            plist->push_back(parton);                             // add the parton to the main list
          }
      }
          eventNumber++;
        }
    }
    }

      //Now, determine the partners of each heavy quark:
      for(int i = 0; i < plist->size(); i++){
  if(i%2 == 0){
    plist->at(i).antiI(i+1);
  }
  else
    plist->at(i).antiI(i-1);
      }

      totalNNs = eventNumber;
      //output for plot
      ofstream fout1("./output/density.dat",ios::out); 
      ofstream fout2("./output/densityB.dat",ios::out); 
      ofstream fout3("./output/colls.dat",ios::out); 
      ofstream fout4("./output/NN.dat",ios::app); 
      
      for (int i = 0; i < ixmax; i++)
  for (int j = 0; j < ixmax; j++)
    {
      fout1 <<  i << " " << j << " " << nucALat[i][j] << endl;
      if ( j == ixmax-1 ) fout1 << endl;
    }
      
      for (int i = 0; i < ixmax; i++)
  for (int j = 0; j < ixmax; j++)
    {
      fout2 <<  i << " " << j << " " << nucBLat[i][j] << endl;
      if ( j == ixmax-1 ) fout2 << endl;
    }
      
      for (int i = 0; i < ixmax; i++)
  for (int j = 0; j < ixmax; j++)
    {
      fout3 <<  i << " " << j << " " << collLat[i][j] << endl;
      if ( j == ixmax-1 ) fout3 << endl;
    }
      
      fout4 << numberOfpp << endl;
      
      fout1.close();
      fout2.close();
      fout3.close();
      fout4.close();
      return eventNumber;
    }
  else if (fullEvent==0) // sample only one hard collision per heavy-ion event
    {
      Parton parton;                                                // new Parton object
      int n1;
      int n2;
      int Projectile_A = glauber->get_Projetile_A();
      int Target_A = glauber->get_Target_A();
      int Projectile_Z = glauber->get_Projetile_Z();
      int Target_Z = glauber->get_Target_Z();
      if ( random->genrand64_real1() > (double)Projectile_Z/(double)Projectile_A ) n1=2112;
      else n1=2212;
      if ( random->genrand64_real1() > (double)Target_Z/(double)Target_A ) n2=2112;
      else n2=2212;

      int pythiaWorked = 0;
      while(pythiaWorked == 0){
  pythia.next();                                                // generate event with pythia
//  pythiaWorked = pythia.next(n1,n2);                                                // generate event with pythia
      }

      ReturnValue rv;
      //if (evolution && fixedTemperature==0) 
      if (!allFromCenter)
  {
    if(nbinFromFile == 1){
      //Sample a random element of the list of number of collisions. -CFY 11/2/2010
      double randomr = random->genrand64_real2();
      int randomint = (int)(randomr*((double)Nbin));
      rv.x = binary[randomint][0];
      rv.y = binary[randomint][1];
    }
    else if (glauberEnvelope)
      rv=glauber->SamplePAB(random);          // determine position in x-y plane from Glauber model with Metropolis
    else
      rv=glauber->SamplePABRejection(random); // determine position in x-y plane from Glauber model with rejection method
    cout << "rv.x = " << rv.x << ", rv.y = " << rv.y << endl;
  }
      else
  {
    rv.x=0.;
    rv.y=0.;
  }
      if (trackHistory)
  {
    rv.x=initialXjet;
    rv.y=initialYjet;
  }
      for (int i = 0; i < pythia.event.size(); ++i) 
  {
    if (pythia.event[i].isFinal() && (abs(pythia.event[i].id()) == 4 || abs(pythia.event[i].id()) == 5) )
      {
        double dx, dy, dz;
        dx = dy = dz = 0.;
        if(abs(pythia.event[i].id() ) == 4){
    dx = gsl_ran_gaussian(gsl_rand, charmWidth);
    dy = gsl_ran_gaussian(gsl_rand, charmWidth);
        }
        if(abs(pythia.event[i].id() ) == 5){
    dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
    dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
        }

        if( moveBeforeTau0 == 1){
    Vec4 partonP = pythia.event[i].p();
    double pzP = partonP.pz();
    double EP = sqrt(pythia.event[i].m()*pythia.event[i].m()
         +partonP.px()*partonP.px()
         +partonP.py()*partonP.py()
         +partonP.pz()*partonP.pz() );
    double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
    dx += partonP.px()*dt/EP;
    dy += partonP.py()*dt/EP;
    dz += partonP.pz()*dt/EP;
        }

        parton.id(pythia.event[i].id());                      // set parton id
        parton.status(pythia.event[i].status());              // set parton status
        parton.mass(pythia.event[i].m());                // set mass
                    
        parton.x(rv.x+dx);                                   // set position
        parton.y(rv.y+dy);
//        if (fixedTemperature==0)
//    parton.tini(hydroTau0);
//        else
//    parton.tini(0.);
//        parton.xini(rv.x+dx);                                // set initial position to remember
//        parton.yini(rv.y+dy);
                    
        parton.z(dz);
//        parton.zini(dz);
        //Set the old positions to be the current ones, for the first timestep:
        parton.xOld(rv.x+dx);                                   // set position
        parton.yOld(rv.y+dy);
        parton.zOld(dz);
        parton.col(pythia.event[i].col());                 // set color 
        parton.acol(pythia.event[i].acol());                  // set anti-color
        parton.frozen(0);                                     // parton is not frozen (will evolve)
        parton.formed(0);                                     // parton is not formed immediately (will evolve)
        parton.inhydro(0);                                     // parton is not formed immediately (will evolve)
        if(setInitialPToZero == 0){
    parton.p(pythia.event[i].p());                        // set momentum
    parton.Initialp(pythia.event[i].p());                        // set momentum
    parton.pOld(pythia.event[i].p());                        // set momentum
//    parton.pini(pythia.event[i].p());                     // set momentum
        }
        else{
    parton.p(0.,0.,0.);                                   // set momentum
    parton.Initialp(0.,0.,0.);                                   // set momentum
    parton.pOld(0.,0.,0.);                                   // set momentum
//    Vec4 pinitial(pythia.event[i].mass(), 0.,0.,0.);
//    parton.pini(pinitial);                               // set momentum
        }
        parton.tFinal(0.);                                    // Set the freezeout time initially to zero
//        parton.elasticCollisions(0);                          // set initial no. of el colls to zero
        //parton.splits(0);                                   // set initial no. of splittings to zero
        parton.source(0);                                     // All intial partons have itsSource=0
        plist->push_back(parton);                             // add the parton to the main list
      }
  }

      //Now, determine the partners of each heavy quark:
      for(int i = 0; i < plist->size(); i++){
  if(i%2 == 0){
    plist->at(i).antiI(i+1);
  }
  else
    plist->at(i).antiI(i-1);
      }

      return 1;
    }
}

int MARTINI::generateEventHeavyQuarks_w_colllist_OR_IPGfile(vector<Parton> *plist)
{
  pythia.event.reset();
  if (fullEvent==1)
  {
    Parton parton;                               // new Parton object
    int n1=0;
    int n2=0;
    string dummy;
    vector<ReturnValue> BinaryColl;

    if(nbinFromFile == 1)
    {
      ifstream fin(binary_filename.c_str());
      if(!fin)
      {
        cout << "ERROR: " << binary_filename << " not found." << endl;
        exit(1);
      }
      while (true) {
        ReturnValue rv;
        double collx, colly;
        int dum;
        getline(fin, dummy);
        stringstream ss(dummy);;

        ss >> collx >> colly;
        if (fin.eof()) break;

        rv.x = collx;
        rv.y = colly;
        BinaryColl.push_back(rv);
      }
    } else
    {
      ifstream fin(IPG_filename.c_str());
      if(!fin)
      {
        cout << "ERROR: " << IPG_filename << " not found." << endl;
        exit(1);
      }
      double **ener_density,**cdf_xy, *cdf_x, dx, dy;
      int nx, ny;

      fin  >> dummy >> dummy >> dummy >> nx >> dummy >> ny
           >> dummy >> dummy >> dummy >> dx >> dummy >> dy;
      ener_density = new double* [nx];
      cdf_xy = new double* [nx];
      cdf_x  = new double[nx];
      for (int ix = 0; ix < nx; ix++)
      {
        ener_density[ix] = new double[ny];
        cdf_xy[ix] = new double[ny];
      }

      double tot_ener_dens = 0.;
      double density;
      for (int ix = 0; ix < nx; ix++) {
       for (int iy = 0; iy < ny; iy++) {
            std::getline(fin, dummy);
            std::stringstream ss(dummy);
            ss >> dummy >> dummy >> dummy
               >> density >> dummy >> dummy >> dummy >> dummy
               >> dummy >> dummy >> dummy >> dummy >> dummy
               >> dummy >> dummy >> dummy >> dummy >> dummy;
            ener_density[ix][iy] = density;
            tot_ener_dens += density;
        }
        std::getline(fin, dummy);
      }
      fin.close();

      for (int ix = 0; ix < nx; ix++) {
         double tot_dens = 0.;
         for (int iy = 0; iy < ny; iy++)
           tot_dens += ener_density[ix][iy];
         cdf_x[ix] = tot_dens;
         cdf_xy[ix][0] = ener_density[ix][0]/(tot_dens+1.e-10);
         for (int iy = 1; iy < ny; iy++)
           cdf_xy[ix][iy] = cdf_xy[ix][iy-1] + ener_density[ix][iy]/(tot_dens+1.e-10);
      }

      cdf_x[0] /= tot_ener_dens;
      for (int ix = 1; ix < nx; ix++)
        cdf_x[ix] = cdf_x[ix-1] + cdf_x[ix]/tot_ener_dens;
      
      for (int ii = 0; ii < NColl; ii++)
      {
        ReturnValue rv;
        double collx, colly;
        double rand1 = random->genrand64_real1();
        double rand2 = random->genrand64_real1();
        for (int ix = 0; ix < nx; ix++)
        {
          int icoll_x;
          if (rand1 >= cdf_x[ix] && rand1 < cdf_x[ix+1] )
          {
            icoll_x = ix;
            int icoll_y;
            for (int iy = 0; iy < ny; iy++)
            {
              if (rand2 >= cdf_xy[icoll_x][iy] && rand2 <= cdf_xy[icoll_x][iy+1])
              {
                icoll_y = iy;
                 collx = (-0.5*((double)nx)+((double)icoll_x))*dx;
                 collx += dx*(rand1 - cdf_x[icoll_x]);
                 colly = (-0.5*((double)ny)+((double)icoll_y))*dy;
                 colly += dx*(rand2 - cdf_xy[icoll_x][icoll_y]);
                 break;
              }
            }
            break;
          }
        }
        rv.x = collx;
        rv.y = colly;
        BinaryColl.push_back(rv);
      }
      for (int ix = 0; ix < nx; ix++)
      {
        delete[] ener_density[ix];
        delete[] cdf_xy[ix];
      }
      delete[] ener_density;
      delete[] cdf_xy;
      delete[] cdf_x;
      //while (true) {
      //  ReturnValue rv;
      //  double collx, colly;
      //  int dum;
      //  getline(fin, dummy);
      //  stringstream ss(dummy);;

      //  ss >> dum >> collx >> colly;
      //  if (fin.eof()) break;

      //  rv.x = collx;
      //  rv.y = colly;
      //  BinaryColl.push_back(rv);
      //}
    }

    int done;
    int numberOfpp;
    int eventNumber=0;
    numberOfpp=0;
    int number_of_binary_collisions = BinaryColl.size();
//    cout << "Number of binary collisions = " << number_of_binary_collisions << endl;
    int Projectile_A=glauber->get_Projetile_A();
    int Target_A=glauber->get_Projetile_A();
    int Projectile_Z = glauber->get_Projetile_Z();
    int Target_Z = glauber->get_Projetile_Z();
    for (int icoll = 0; icoll < number_of_binary_collisions; icoll++)
    { 
      if ( random->genrand64_real1() < totalHQXSec/inelasticXSec )
      {
        if ( random->genrand64_real1() > (double)Projectile_Z/(double)Projectile_A ) n1=2112;
        else n1=2212;
        if ( random->genrand64_real1() > (double)Target_Z/(double)Target_A ) n2=2112;
        else n2=2212;

	int pythiaWorked = 0;
        while(pythiaWorked == 0){
          pythia.next();
//          pythiaWorked = pythia.next(n1,n2);
        }                                                // generate event with pythia

        done = 0;
        for (int ip = 0; ip < pythia.event.size(); ++ip) 
        {
          if (pythia.event[ip].status()>0 && (abs(pythia.event[ip].id())==4 || abs(pythia.event[ip].id()) == 5) )
          // If the parton is final, and is a heavy quark, then put it into the list.
          {
            parton.id(pythia.event[ip].id());                     // set parton id
            parton.status(pythia.event[ip].status());             // set parton status
            parton.mass(pythia.event[ip].m());                // set mass
            if (fixedTemperature==0)
            {
              //Random numbers added to the initial positions:
              double dx, dy, dz;
              dx = dy = dz = 0.;
              if(abs(parton.id()) == 4 && charmWidth > 0.01)
              {
                dx = gsl_ran_gaussian(gsl_rand, charmWidth);
                dy = gsl_ran_gaussian(gsl_rand, charmWidth);
              }
              if(abs(parton.id()) == 5 && bottomWidth > 0.01)
              {
                dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
                dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
              }

              if( moveBeforeTau0 == 1)
              {
                Vec4 partonP = pythia.event[ip].p();
                double pzP = partonP.pz();
                double EP = sqrt(pythia.event[ip].m()*pythia.event[ip].m()
                     +partonP.px()*partonP.px()
                     +partonP.py()*partonP.py()
                     +partonP.pz()*partonP.pz() );
                double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
                dx += partonP.px()*dt/EP;
                dy += partonP.py()*dt/EP;
                dz += partonP.pz()*dt/EP;
              }

              parton.x(BinaryColl[icoll].x+dx);       // set position
              parton.y(BinaryColl[icoll].y+dy);
              parton.z(dz);
	      parton.prevtau(0.);
              //Set the old positions to be the current ones, for the first timestep:
              parton.xOld(BinaryColl[icoll].x+dx);       // set position
              parton.yOld(BinaryColl[icoll].y+dy);
              parton.zOld(dz);
              if ( done == 0 )
              {
                numberOfpp++;
                done = 1;
              }
            }
            else
            {
              double dx, dy, dz;
              dx = dy = dz = 0.;
              if(abs(parton.id()) == 4)
              {
                dx = gsl_ran_gaussian(gsl_rand, charmWidth);
                dy = gsl_ran_gaussian(gsl_rand, charmWidth);
              }
              if(abs(parton.id()) == 5)
              {
                dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
                dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
              }

              if( moveBeforeTau0 == 1)
              {
                Vec4 partonP = pythia.event[ip].p();
                double pzP = partonP.pz();
                double EP = sqrt(pythia.event[ip].m()*pythia.event[ip].m()
                     +partonP.px()*partonP.px()
                     +partonP.py()*partonP.py()
                     +partonP.pz()*partonP.pz() );
                double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
                dx += partonP.px()*dt/EP;
                dy += partonP.py()*dt/EP;
                dz += partonP.pz()*dt/EP;
              }

              parton.x(dx);     
              parton.y(dy);
              parton.z(dz);
              //Set the old positions to be the current ones, for the first timestep:
              parton.xOld(dx);     
              parton.yOld(dy);
              parton.zOld(dz);
            }
            parton.col(pythia.event[ip].col());                 // set color 
            parton.acol(pythia.event[ip].acol());                 // set anti-color
            parton.frozen(0);                                     // parton is not frozen (will evolve)
            parton.formed(0);                                     // parton is not immediately formed
            parton.inhydro(0);                                     // parton is not immediately formed
            if(setInitialPToZero == 0)
            {
              parton.p(pythia.event[ip].p());  // set momentum
              parton.Initialp(pythia.event[ip].p());  // set momentum
              parton.pOld(pythia.event[ip].p());  // set momentum
            }
            else
            {
              parton.p(0.,0.,0.);
              parton.Initialp(0.,0.,0.);
              parton.pOld(0.,0.,0.);
            }
            parton.tFinal(0.);                                    // Set the initial freezeout time to zero
            parton.source(0);                                     // All initial partons have itsSource=0
            plist->push_back(parton);                             // add the parton to the main list
          }
        }
        eventNumber++;
      }
    }

    //Now, determine the partners of each heavy quark:
    for(int i = 0; i < plist->size(); i++)
    {
      if(i%2 == 0)
      {
        plist->at(i).antiI(i+1);
      }
      else
      plist->at(i).antiI(i-1);
    }

    totalNNs = eventNumber;
    ofstream fout4("./output/NN.dat",ios::app); 
      
    fout4 << numberOfpp << endl;
      
    fout4.close();
    return eventNumber;
  }
  else // sample only one hard collision per heavy-ion event
  {
    cout << "ERROR: generateEventHeavyQuarks_withBinaryCollisionsList works only for fill events and not for only one HQ per collision" << endl;
    exit(1);
  }
}
//int MARTINI::generateEventHeavyQuarks(vector<Parton> *plist)
//{
//  pythia.event.reset();
//
//  if (fullEvent==1)
//    {
//      // Sample Nu (which gives the number of nucleons at radial position r) for both nuclei
//      // then lay them on top of each other and sample the number of charm and bottom events with 
//      // area = totalHQXSec.
//      
//      //cout << "Number of binary collisions = " << glauber->TAB() << endl; 
//      //cout << " A=" << glauber->nucleusA() << endl;
//
//      double A;
//      double Z;
//      Parton parton;                               // new Parton object
//      double b = glauberImpactParam;
//      int posx;
//      int posy;
//      int n1=0;
//      int n2=0;
////      double r = 1.2*pow(glauber->nucleusA(),1./3.);
//      double r = 1.2*pow(nucleusA(),1./3.);
//
//      nucleusA.clear();
//      nucleusB.clear();
//      sampleTA();                                  // populate the lists nucleusA and nucleusB with position data of the nucleons
//
//      double cellLength = sqrt(inelasticXSec*0.1); // compute cell length in fm (1 mb = 0.1 fm^2)
//      double latticeLength = 4.*r;                 // spread the lattice 2r in both (+/-) directions (total length = 4r)
//      int ixmax = ceil(latticeLength/cellLength);
//      ixmax*=2;
//      double xmin = -latticeLength;
//      
//      int nucALat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus A
//      int nucBLat[ixmax][ixmax];                   // lattice with ixmax*ixmax cells for nucleus B
//      int collLat[ixmax][ixmax];                   // lattice with the positions of the interactions -> for information only
//
//      int countCollisionsA[100];
//      int countCollisionsB[100];
//
//      for (int i=0; i<100; i++)
//  {
//    countCollisionsA[i]=0;
//    countCollisionsB[i]=0;
//  }
//
//      for (int i = 0; i < ixmax; i++)              // initialize cells: zero nucleons in every cell
//  for (int j = 0; j < ixmax; j++)
//    {
//      nucALat[i][j] = 0;
//      nucBLat[i][j] = 0;
//      collLat[i][j] = 0;
//    }
//
//      for (int i = 0; i < nucleusA.size(); i++)    // fill nucleons into cells
//  {
//    posx = floor((nucleusA.at(i).x-xmin-b/2.)/cellLength);
//    posy = floor((nucleusA.at(i).y-xmin)/cellLength);
//    nucALat[posx][posy] += 1;
//        
//    posx = floor((nucleusB.at(i).x-xmin+b/2.)/cellLength);
//    posy = floor((nucleusB.at(i).y-xmin)/cellLength);
//    nucBLat[posx][posy] += 1;
//
//  }
//
//      int done;
//      int numberOfpp;
//      int eventNumber=0;
//      numberOfpp=0;
//      //cout << "jetXSec=" << jetXSec << endl;
//      //cout << "inelasticXSec=" << inelasticXSec << endl;
//      //A=glauber->nucleusA();
//      A=nucleusA();
//      if (A==197.) //gold - Au 
//  Z=79.; 
//      else if (A==63.) //copper - Cu
//  Z=29.;
//      else if (A==208.) //lead - Pb
//  Z=82.;
//      else
//  {
//    Z=0;
//  }
//      for (int ix = 0; ix < ixmax; ix++) 
//  for (int iy = 0; iy < ixmax; iy++)
//    {
//      for (int i = 0; i < nucALat[ix][iy]; i++)
//        for (int j = 0; j < nucBLat[ix][iy]; j++)
//    {
//      //cout << "ix+b/2=" << static_cast<int>(ix+b/2.) << " # there=" << nucALat[static_cast<int>(ix+b/2.)][iy] << endl;
//      if ( random->genrand64_real1() < totalHQXSec/inelasticXSec )
//        {
//          countCollisionsA[i]+=1;
//          countCollisionsB[j]+=1;
//          double xPositionInCell = (random->genrand64_real1())*cellLength;
//          double yPositionInCell = (random->genrand64_real1())*cellLength;
//          // see if colliding nucleons are neutrons
//          if ( random->genrand64_real1() > Z/A ) n1=1;
//          else n1=0;
//          if ( random->genrand64_real1() > Z/A ) n2=1;
//          else n2=0;
//
//          int pythiaWorked = 0;
//          while(pythiaWorked == 0){
//      pythiaWorked = pythia.next(n1,n2);
//          }                                                // generate event with pythia
//
//          done = 0;
//          for (int ip = 0; ip < pythia.event.size(); ++ip) 
//      {
//        if (pythia.event[ip].status()>0 && (abs(pythia.event[ip].id())==4 || abs(pythia.event[ip].id()) == 5) )
//          // If the parton is final, and is a heavy quark, then put it into the list.
//          {
//            parton.id(pythia.event[ip].id());                     // set parton id
//            parton.status(pythia.event[ip].status());             // set parton status
//            parton.mass(pythia.event[ip].m());                // set mass
//            if (fixedTemperature==0)
//        {
//          //Random numbers added to the initial positions:
//          double dx, dy, dz;
//          dx = dy = dz = 0.;
//          if(abs(parton.id()) == 4 && charmWidth > 0.01){
//            dx = gsl_ran_gaussian(gsl_rand, charmWidth);
//            dy = gsl_ran_gaussian(gsl_rand, charmWidth);
//          }
//          if(abs(parton.id()) == 5 && bottomWidth > 0.01){
//            dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
//            dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
//          }
//
//          if( moveBeforeTau0 == 1){
//            Vec4 partonP = pythia.event[ip].p();
//            double pzP = partonP.pz();
//            double EP = sqrt(pythia.event[ip].m()*pythia.event[ip].m()
//                 +partonP.px()*partonP.px()
//                 +partonP.py()*partonP.py()
//                 +partonP.pz()*partonP.pz() );
//            double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
//            dx += partonP.px()*dt/EP;
//            dy += partonP.py()*dt/EP;
//            dz += partonP.pz()*dt/EP;
//          }
//
//          parton.x(xmin+(ix)*cellLength+xPositionInCell+dx);       // set position
//          parton.y(xmin+(iy)*cellLength+yPositionInCell+dy);
//          parton.xini(xmin+(ix)*cellLength+xPositionInCell+dx);       // set position
//          parton.yini(xmin+(iy)*cellLength+yPositionInCell+dy);
//          parton.z(dz);
//          parton.zini(dz);
//          parton.tini(hydroTau0);
//          //Set the old positions to be the current ones, for the first timestep:
//          parton.xOld(xmin+(ix)*cellLength+xPositionInCell+dx);       // set position
//          parton.yOld(xmin+(iy)*cellLength+yPositionInCell+dy);
//          parton.zOld(dz);
//          if ( done == 0 )
//            {
//              posx = floor((parton.x()-xmin)/cellLength);
//              posy = floor((parton.y()-xmin)/cellLength);
//              collLat[posx][posy] += 1;
//                          
//              //cout << numberOfpp << endl;
//              numberOfpp++;
//              done = 1;
//            }
//        }
//            else
//        {
//          double dx, dy, dz;
//          dx = dy = dz = 0.;
//          if(abs(parton.id()) == 4){
//            dx = gsl_ran_gaussian(gsl_rand, charmWidth);
//            dy = gsl_ran_gaussian(gsl_rand, charmWidth);
//          }
//          if(abs(parton.id()) == 5){
//            dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
//            dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
//          }
//
//          if( moveBeforeTau0 == 1){
//            Vec4 partonP = pythia.event[ip].p();
//            double pzP = partonP.pz();
//            double EP = sqrt(pythia.event[ip].m()*pythia.event[ip].m()
//                 +partonP.px()*partonP.px()
//                 +partonP.py()*partonP.py()
//                 +partonP.pz()*partonP.pz() );
//            double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
//            dx += partonP.px()*dt/EP;
//            dy += partonP.py()*dt/EP;
//            dz += partonP.pz()*dt/EP;
//          }
//
//          parton.x(dx);     
//          parton.y(dy);
//          parton.xini(dx);     
//          parton.yini(dy);
//          parton.z(dz);
//          parton.zini(dz);
//          parton.tini(0.);
//          //Set the old positions to be the current ones, for the first timestep:
//          parton.xOld(dx);     
//          parton.yOld(dy);
//          parton.zOld(dz);
//        }
//            parton.col(pythia.event[ip].col());                 // set color 
//            parton.acol(pythia.event[ip].acol());                 // set anti-color
//            parton.frozen(0);                                     // parton is not frozen (will evolve)
//            if(setInitialPToZero == 0){
//        parton.p(pythia.event[ip].p());  // set momentum
//        parton.pOld(pythia.event[ip].p());  // set momentum
//        parton.pini(pythia.event[ip].p());                    // set initial momentum
//            }
//            else{
//        parton.p(0.,0.,0.);
//        parton.pOld(0.,0.,0.);
//        Vec4 pinitial(pythia.event[ip].mass(),0.,0.,0.);
//        parton.pini(pinitial);                    // set initial momentum
//            }
//            parton.tFinal(0.);                                    // Set the initial freezeout time to zero
//            parton.eventNumber(eventNumber);
//            parton.elasticCollisions(0);                          // set initial no. of el colls to zero
//            parton.source(0);                                     // All initial partons have itsSource=0
//            plist->push_back(parton);                             // add the parton to the main list
//          }
//      }
//          eventNumber++;
//        }
//    }
//    }
//
//      //Now, determine the partners of each heavy quark:
//      for(int i = 0; i < plist->size(); i++){
//  if(i%2 == 0){
//    plist->at(i).antiI(i+1);
//  }
//  else
//    plist->at(i).antiI(i-1);
//      }
//
//      totalNNs = eventNumber;
//      //output for plot
//      ofstream fout1("./output/density.dat",ios::out); 
//      ofstream fout2("./output/densityB.dat",ios::out); 
//      ofstream fout3("./output/colls.dat",ios::out); 
//      ofstream fout4("./output/NN.dat",ios::app); 
//      
//      for (int i = 0; i < ixmax; i++)
//  for (int j = 0; j < ixmax; j++)
//    {
//      fout1 <<  i << " " << j << " " << nucALat[i][j] << endl;
//      if ( j == ixmax-1 ) fout1 << endl;
//    }
//      
//      for (int i = 0; i < ixmax; i++)
//  for (int j = 0; j < ixmax; j++)
//    {
//      fout2 <<  i << " " << j << " " << nucBLat[i][j] << endl;
//      if ( j == ixmax-1 ) fout2 << endl;
//    }
//      
//      for (int i = 0; i < ixmax; i++)
//  for (int j = 0; j < ixmax; j++)
//    {
//      fout3 <<  i << " " << j << " " << collLat[i][j] << endl;
//      if ( j == ixmax-1 ) fout3 << endl;
//    }
//      
//      fout4 << numberOfpp << endl;
//      
//      fout1.close();
//      fout2.close();
//      fout3.close();
//      fout4.close();
//      return eventNumber;
//    }
//  else if (fullEvent==0) // sample only one hard collision per heavy-ion event
//    {
//      Parton parton;                                                // new Parton object
//      double A;
//      double Z;
//      int n1;
//      int n2;
//      //A=glauber->nucleusA();
//      A=nucleusA();
//      if (A==197.) //gold - Au 
//  Z=79.; 
//      else if (A==63.) //copper - Cu
//  Z=29.;
//      else if (A==208.) //lead - Pb
//  Z=82.;
//      else
//  {
//    Z=0;
//  }
//      if ( random->genrand64_real1() > Z/A ) n1=1;
//      else n1=0;
//      if ( random->genrand64_real1() > Z/A ) n2=1;
//      else n2=0;
//
//      int pythiaWorked = 0;
//      while(pythiaWorked == 0){
//  pythiaWorked = pythia.next(n1,n2);                                                // generate event with pythia
//      }
//
//      ReturnValue rv;
//      //if (evolution && fixedTemperature==0) 
//      if (!allFromCenter)
//  {
//    if(nbinFromFile == 1){
//      //Sample a random element of the list of number of collisions. -CFY 11/2/2010
//      double randomr = random->genrand64_real2();
//      int randomint = (int)(randomr*((double)Nbin));
//      rv.x = binary[randomint][0];
//      rv.y = binary[randomint][1];
//    }
//    else if (glauberEnvelope)
//      rv=glauber->SamplePAB(random);          // determine position in x-y plane from Glauber model with Metropolis
//    else
//      rv=glauber->SamplePABRejection(random); // determine position in x-y plane from Glauber model with rejection method
//    cout << "rv.x = " << rv.x << ", rv.y = " << rv.y << endl;
//  }
//      else
//  {
//    rv.x=0.;
//    rv.y=0.;
//  }
//      if (trackHistory)
//  {
//    rv.x=initialXjet;
//    rv.y=initialYjet;
//  }
//      for (int i = 0; i < pythia.event.size(); ++i) 
//  {
//    if (pythia.event[i].isFinal() && (abs(pythia.event[i].id()) == 4 || abs(pythia.event[i].id()) == 5) )
//      {
//        double dx, dy, dz;
//        dx = dy = dz = 0.;
//        if(abs(pythia.event[i].id() ) == 4){
//    dx = gsl_ran_gaussian(gsl_rand, charmWidth);
//    dy = gsl_ran_gaussian(gsl_rand, charmWidth);
//        }
//        if(abs(pythia.event[i].id() ) == 5){
//    dx = gsl_ran_gaussian(gsl_rand, bottomWidth);
//    dy = gsl_ran_gaussian(gsl_rand, bottomWidth);
//        }
//
//        if( moveBeforeTau0 == 1){
//    Vec4 partonP = pythia.event[i].p();
//    double pzP = partonP.pz();
//    double EP = sqrt(pythia.event[i].m()*pythia.event[i].m()
//         +partonP.px()*partonP.px()
//         +partonP.py()*partonP.py()
//         +partonP.pz()*partonP.pz() );
//    double dt = hydroTau0/sqrt(1.-(pzP/EP)*(pzP/EP));
//    dx += partonP.px()*dt/EP;
//    dy += partonP.py()*dt/EP;
//    dz += partonP.pz()*dt/EP;
//        }
//
//        parton.id(pythia.event[i].id());                      // set parton id
//        parton.status(pythia.event[i].status());              // set parton status
//        parton.mass(pythia.event[i].m());                // set mass
//                    
//        parton.x(rv.x+dx);                                   // set position
//        parton.y(rv.y+dy);
//        if (fixedTemperature==0)
//    parton.tini(hydroTau0);
//        else
//    parton.tini(0.);
//        parton.xini(rv.x+dx);                                // set initial position to remember
//        parton.yini(rv.y+dy);
//                    
//        parton.z(dz);
//        parton.zini(dz);
//        //Set the old positions to be the current ones, for the first timestep:
//        parton.xOld(rv.x+dx);                                   // set position
//        parton.yOld(rv.y+dy);
//        parton.zOld(dz);
//        parton.col(pythia.event[i].col());                 // set color 
//        parton.acol(pythia.event[i].acol());                  // set anti-color
//        parton.frozen(0);                                     // parton is not frozen (will evolve)
//        if(setInitialPToZero == 0){
//    parton.p(pythia.event[i].p());                        // set momentum
//    parton.pOld(pythia.event[i].p());                        // set momentum
//    parton.pini(pythia.event[i].p());                     // set momentum
//        }
//        else{
//    parton.p(0.,0.,0.);                                   // set momentum
//    parton.pOld(0.,0.,0.);                                   // set momentum
//    Vec4 pinitial(pythia.event[i].mass(), 0.,0.,0.);
//    parton.pini(pinitial);                               // set momentum
//        }
//        parton.tFinal(0.);                                    // Set the freezeout time initially to zero
//        parton.elasticCollisions(0);                          // set initial no. of el colls to zero
//        //parton.splits(0);                                   // set initial no. of splittings to zero
//        parton.source(0);                                     // All intial partons have itsSource=0
//        plist->push_back(parton);                             // add the parton to the main list
//      }
//  }
//
//      //Now, determine the partners of each heavy quark:
//      for(int i = 0; i < plist->size(); i++){
//  if(i%2 == 0){
//    plist->at(i).antiI(i+1);
//  }
//  else
//    plist->at(i).antiI(i-1);
//      }
//
//      return 1;
//    }
//}

// Evolves heavy quarks according to Langevin dynamics. -CFY 4/10/2011/
int MARTINI::evolveHeavyQuarks(vector<Parton>  **plist, int counter, int it)
{
  
  cout.precision(6);
  
  HydroInfo hydroInfo, antihydroInfo;
  HydroInfoViscous hydroViscousInfo, antihydroViscousInfo;
  
  Vec4 vecp, antivecp;                          // momentum four-vectors
  double pmu[4], pmuRest[4], pmuRestNew[4], pmuNew[4]; //Easy to Lorentz-boost arrays with my subroutine 
  double antipmu[4], antipmuRest[4], antipmuRestNew[4], antipmuNew[4]; //Easy to Lorentz-boost arrays with my subroutine 
  
  int imax = plist[0]->size();        // number of partons in the list
  int id, antiI; //newID;                      // original and new parton's ID
  //int col, acol, newCol, newAcol;     // number of the color string attached to the old, and new parton
  int ix, iy, iz, itau;               // parton's position in cell coordinates
  int antiix, antiiy, antiiz, antiitau; // parton's position in cell coordinates
  int ixmax, izmax;                   // maximum cell number in x- and y-direction, and z-direction
  
  double dt = dtfm/hbarc;             // time step in [GeV^(-1)] - hbarc is defined in Basics.h 
  // Gamma = \int dGamma dp is given in units of GeV, so dt*Gamma is dimensionless
  double dtflow, antidtflow;                      // dt in the fluid's rest frame
  double T, antiT, T_ave;             // temperature
  double M, antiM;                    // Parton's mass [GeV]
  double x, y, z;                     // parton's position in [fm]
  double antix, antiy, antiz       ;  // parton's position in [fm]
  double t, tau, antitau;             // lab time and lab tau
  double vx, vy, vz;                  // flow velocity of the cell
  double antivx, antivy, antivz;      // flow velocity of the cell
  double beta;                        // absolute value of the flow velocity
  double antibeta;                        // absolute value of the flow velocity
  double gamma;                       // gamma factor 
  double antigamma;                       // gamma factor 
  double umu[4], umub[4];             // Fluid's 4-velocity
  double antiumu[4], antiumub[4];             // Fluid's 4-velocity
  const double rd=12.3;               // ratio of degrees of freedom g_QGP/g_H (for Kolb hydro only)
  double Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta;
  double Wyy, Wyeta, Wetaeta, BulkPi, Entropy, cs2;
  
  double kappa_T, kappa_L, sigma, eta;
  double kicktriplet[3];
  double r, force, F_x, F_y, F_z;
  
  ixmax = static_cast<int>(2.*hydroXmax/hydroDx+0.0001);
  
  if ( hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6) izmax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);
  
  //cout << "counter = " << counter << ", it = " << it << endl;
  
  for ( int i=0; i<imax; i++)                                   // loop over all partons
    { 
      //cout << "plist[0]->at(" << i << ").frozen() = " << plist[0]->at(i).frozen() << endl;
      if(plist[0]->at(i).frozen() == 0){
  counter++;                                                // counter that characterizes every step uniquely
  id = plist[0]->at(i).id();                                // id of parton i 
  
  if(abs(id)!=4 && abs(id)!=5) continue;                    // We immediately leave the loop of the parton is not a heavy quark!
  
  
  t = it*dtfm;        //had a +tau0 here. now let partons travel doing their vacuum shower stuff first...
  
  vecp = plist[0]->at(i).p();                               // four-vector p of parton i
  M = plist[0]->at(i).mass();                               // Mass of the parton (only used for charm and bottom quarks)
  pmu[1] = vecp.px();
  pmu[2] = vecp.py();
  pmu[3] = vecp.pz();
  pmu[0] = sqrt(M*M+pmu[1]*pmu[1]+pmu[2]*pmu[2]+pmu[3]*pmu[3]); //vecp loaded into an array for ease with LorentzBooster
  
  x = plist[0]->at(i).x();                                  // x value of position in [fm]
  y = plist[0]->at(i).y();                                  // y value of position in [fm]
  z = plist[0]->at(i).z();                                  // z value of position in [fm]
  
  // boost - using original position...
  ix = floor((hydroXmax+x)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
  iy = floor((hydroXmax+y)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
  // note that x and y run from -hydroXmax to +hydroXmax
  // and ix and iy from 0 to 2*hydroXmax
  // hence the (hydroXmax+x or y) for both
  if (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) 
    {
      iz = floor((hydroZmax+z)/hydroDz+0.0001);             // z-coordinate of the cell we are in now
    }
  
  if (fixedTemperature == 0)
    {
      tau = 0.;           //just give tau some value - will change below
      
      if (z*z<t*t) tau = sqrt(t*t-z*z);                     // determine tau (z*z) should always be less than (t*t)
      //cout << "tau = " << tau << endl;
      if (tau<hydroTau0||tau>hydroTauMax-1) continue;       // do not evolve if tau<tau0 or tau>tauMax
      
      // get the temperature and flow velocity at the current position and time:
      hydroInfo = hydroSetup->getHydroValues(x, y, z, t);// hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
               //hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
      
      T = hydroInfo.T;
      //cout << "T = " << T << "for particle " << i << endl;
      vx = hydroInfo.vx;
      vy = hydroInfo.vy;
      vz = hydroInfo.vz;
      
      if(evolutionshear == 1 || evolutionbulk == 1)
            {
              hydroViscousInfo = hydroSetup->getHydroViscousValues(x, y, z, t);
              Wtautau = hydroViscousInfo.Wtautau; 
              Wtaux   = hydroViscousInfo.Wtaux  ; 
              Wtauy   = hydroViscousInfo.Wtauy  ; 
              Wtaueta = hydroViscousInfo.Wtaueta; 
              Wxx     = hydroViscousInfo.Wxx    ; 
              Wxy     = hydroViscousInfo.Wxy    ; 
              Wxeta   = hydroViscousInfo.Wxeta  ; 
              Wyy     = hydroViscousInfo.Wyy    ; 
              Wyeta   = hydroViscousInfo.Wyeta  ; 
              Wetaeta = hydroViscousInfo.Wetaeta; 
              BulkPi  = hydroViscousInfo.BulkPi ; 
              Entropy = hydroViscousInfo.Entropy; 
              cs2     = hydroViscousInfo.cs2    ; 
            }

      //cout << "At t = " << t << " and z = " << z << ", vz = " << vz << endl;
      //cout << "For the charm quark, vz = " << pmu[3]/pmu[0] << endl;
      if (!trackHistory)
        {
    // warn if out of range and move on to next parton (should not happen!)
    if ( ix < 0 || ix >= ixmax ) 
      {
        cout << "WARNING - x out of range" << endl;
        plist[0]->at(i).frozen(1);
        continue;
      }
    if ( iy < 0 || iy >= ixmax )
      {
        cout << "WARNING - y out of range" << endl;
        plist[0]->at(i).frozen(1);
        continue;
      }
    if ( (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) && ( iz < 0 || iz >= izmax ) ) 
      {
        cout << "WARNING - z out of range" << endl;
        plist[0]->at(i).frozen(1);
        continue;
      }
        }
      
      if (T<hydroTfinal && hydroWhichHydro == 2)                 // for 2D hydro evolutions stop at T_final (can change that)
        {
    plist[0]->at(i).tFinal(t);
    continue;                                             // do not evolve if T<T_c
        }
      
      beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of flow velocity in units of c
      gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
      umu[0] = gamma;
      umu[1] = gamma*vx;
      umu[2] = gamma*vy;
      umu[3] = gamma*vz;
      LorentzBooster(umu, pmu, pmuRest);                      // Boost into the fluid's rest frame
      dtflow = dt*pmuRest[0]/pmu[0];                          // dtfm in the fluid's rest frame
      
    }
  else // the fixed temperature case:
    {
      beta = 0.;
      gamma = 1.; 
      T = fixedT;
      umu[0] = 1.; umu[1] = umu[2] = umu[3] = 0.;
      LorentzBooster(umu, pmu, pmuRest);
      dtflow = dt*pmuRest[0]/pmu[0];
    }  
  
  //After all of the various hydro switches, we need to evolve the heavy quarks. This is optimized either for 
  //quarkonium or for single heavy quarks. We differentiate between the two with the switch examineHQ. First, for single 
  //heavy quarks:
  if(examineHQ == 1){
    //Finally evolve the heavy quarks' positions:
    plist[0]->at(i).x(x+vecp.px()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
    plist[0]->at(i).y(y+vecp.py()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
    plist[0]->at(i).z(z+vecp.pz()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
    
    // Sample a triplet of kicks for the heavy quark. For now, we 
    // take the simplest case: constant 2piT*D_c, however later 
    // both kicktriplet[3] and eta will be determined from HTL 
    // results at LO or NLO order:
    if(T > T_C_HQ){
      plist[0]->at(i).tFinal(t); // update the time - finally will be the final time
      double TE = T*sqrt(M*M + pmuRest[1]*pmuRest[1] + pmuRest[2]*pmuRest[2] + pmuRest[3]*pmuRest[3]);
      double pHQ = sqrt(pmuRest[1]*pmuRest[1] + pmuRest[2]*pmuRest[2] + pmuRest[3]*pmuRest[3]);
      if(HQ_rate_type == 0)
      {
      	kappa_T = 4.*M_PI*T*T*T/TwoPiTD_HQ;
        eta = 2.*M_PI*T*T/(TwoPiTD_HQ*M);
        kappa_L = kappa_T;
      } else if (HQ_rate_type == 1)
      {
        double A0, B0, B1;
        hqRates->getTandPdependent_coeff(T, 0., A0, B0, B1);
        eta = A0;
      	kappa_T = 2.*B0;
        kappa_L = kappa_T;
      	//kappa = 2.*eta*T*M;
      } else if (HQ_rate_type == 2)
      {
        double A0, B0, B1;
        hqRates->getTandPdependent_coeff(T, pHQ, A0, B0, B1);
      	kappa_T = 2.*B0;
        eta = A0;
        kappa_L = kappa_T;
      } else if (HQ_rate_type == 3)
      {
        double temp_2piTD;
	if (tau < hydroTau0)
	{
//          temp_2piTD = 20.30405*T + 0.35256;
          temp_2piTD = (10.87937*T - 0.79860)*HQ_Lattice_rate_factor;
	} else
	{
          temp_2piTD = (10.87937*T - 0.79860)*HQ_Lattice_rate_factor;
	}
      	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
      	//kappa_T = 4.*4.*M_PI*T*T*T/temp_2piTD;
        eta = 2.*M_PI*T*T/(temp_2piTD*M);
        kappa_L = kappa_T;

      } else if (HQ_rate_type == 4)
      {
        double temp_2piTD;
	if (tau < hydroTau0)
	{
//          temp_2piTD = 20.30405*T + 0.35256;
          temp_2piTD = (10.87937*T - 0.79860);
	} else
	{
          temp_2piTD = 10.87937*T - 0.79860;
	}
      	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
        eta = 2.*M_PI*T*T/(temp_2piTD*M);
	eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
	kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
        kappa_L = kappa_T;

      } else if (HQ_rate_type == 5)
      {
        double temp_2piTD;
	if (tau < hydroTau0)
	{
          temp_2piTD = (6.07614*T - 0.19409);
	} else
	{
          temp_2piTD = 6.07614*T - 0.19409;
	}
      	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
        eta = 2.*M_PI*T*T/(temp_2piTD*M);
	eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
	kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
        kappa_L = kappa_T;

      } else if (HQ_rate_type == 6)
      {
        double temp_2piTD;
	if (tau < hydroTau0)
	{
          temp_2piTD = (15.682606*T - 1.40311);
	} else
	{
          temp_2piTD = 15.682606*T - 1.40311;
	}
      	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
        eta = 2.*M_PI*T*T/(temp_2piTD*M);
	eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
	kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
        kappa_L = kappa_T;

      } else
      {
      	cout << "Invalid HQ_rate_type" << endl;
	exit(1);
      }

      double kicks[3];
      sigma = sqrt(kappa_T*dtflow);
      kicks[0] = gsl_ran_gaussian(gsl_rand, sigma);
      kicks[1] = gsl_ran_gaussian(gsl_rand, sigma);
      sigma = sqrt(kappa_L*dtflow);
      kicks[2] = gsl_ran_gaussian(gsl_rand, sigma);

      double transv1[3], transv2[3];
      if (pmuRest[3] != 0.)
      {
         transv1[0] = 1.;
         transv1[1] = 1.;
         transv1[2] = -1.*(pmuRest[1]+pmuRest[2])/pmuRest[3];
         double normalization = sqrt(transv1[0]*transv1[0] + transv1[1]*transv1[1] + transv1[2]*transv1[2]);
         for (int ii = 0; ii < 3; ii++) transv1[ii] /= normalization;
      } else
      {
         transv1[0] = 0.;
         transv1[1] = 0.;
         transv1[2] = 1.;
      }

      transv2[0] = pmuRest[2]*transv1[2]/pHQ - pmuRest[3]*transv1[1]/pHQ;
      transv2[1] = pmuRest[3]*transv1[0]/pHQ - pmuRest[1]*transv1[2]/pHQ;
      transv2[2] = pmuRest[1]*transv1[1]/pHQ - pmuRest[2]*transv1[0]/pHQ;

      kicktriplet[0] = kicks[0]*transv1[0] + kicks[1]*transv2[0] + kicks[2]*pmuRest[1]/pHQ;
      kicktriplet[1] = kicks[0]*transv1[1] + kicks[1]*transv2[1] + kicks[2]*pmuRest[2]/pHQ;
      kicktriplet[2] = kicks[0]*transv1[2] + kicks[1]*transv2[2] + kicks[2]*pmuRest[3]/pHQ;

//      kicktriplet[0] = gsl_ran_gaussian(gsl_rand, sigma);
//      kicktriplet[1] = gsl_ran_gaussian(gsl_rand, sigma);
//      kicktriplet[2] = gsl_ran_gaussian(gsl_rand, sigma);
      
      //Update the momentum in the rest frame:
      pmuRestNew[1] = pmuRest[1]*(1.-eta*dtflow)+kicktriplet[0];
      pmuRestNew[2] = pmuRest[2]*(1.-eta*dtflow)+kicktriplet[1];
      pmuRestNew[3] = pmuRest[3]*(1.-eta*dtflow)+kicktriplet[2];
      pmuRestNew[0] = sqrt(M*M+pmuRestNew[1]*pmuRestNew[1]
         +pmuRestNew[2]*pmuRestNew[2]+pmuRestNew[3]*pmuRestNew[3]);
      
      //The flow velocity for boosting back:
      umub[0] = umu[0];
      umub[1] = -umu[1];
      umub[2] = -umu[2];
      umub[3] = -umu[3];
      
      //The boost back into the lab frame:
      LorentzBooster(umub, pmuRestNew, pmuNew);
      
      //Update the particle's momentum:
      plist[0]->at(i).p(pmuNew[1], pmuNew[2], pmuNew[3]);
    }
    else{ plist[0]->at(i).frozen(1);}
      }
  
  //Next, we examine quarkonium:
  if(examineHQ == 2){
    if(id == 4 || id == 5){ //Only the quarks, not the anti-quarks!
      antiI = plist[0]->at(i).antiI();
      //We must duplicate the steps at the beginning of the program for the partner anti-quarks position, momentum, and 
      //other properties:
      antivecp = plist[0]->at(antiI).p();                               // four-vector p of parton i
      antiM = plist[0]->at(antiI).mass();                               // Mass of the parton (only used for charm and bottom quarks)
      antipmu[1] = antivecp.px();
      antipmu[2] = antivecp.py();
      antipmu[3] = antivecp.pz();
      antipmu[0] = sqrt(antiM*antiM+antipmu[1]*antipmu[1]+antipmu[2]*antipmu[2]+antipmu[3]*antipmu[3]); //vecp loaded into an array for ease with LorentzBooster
      
      antix = plist[0]->at(antiI).x();                                  // x value of position in [fm]
      antiy = plist[0]->at(antiI).y();                                  // y value of position in [fm]
      antiz = plist[0]->at(antiI).z();                                  // z value of position in [fm]
      
      antiix = floor((hydroXmax+antix)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
      antiiy = floor((hydroXmax+antiy)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
      if (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) 
        {
    antiiz = floor((hydroZmax+antiz)/hydroDz+0.0001);             // z-coordinate of the cell we are in now
        }
      
      // get the temperature and flow velocity at the current position and time:
      
      if(fixedTemperature == 0){
        antitau = 0.;                                                 // Initially set antitau = 0.
        
        //cout << "antiz = " << antiz << ", t = " << t << endl;
        if (antiz*antiz<t*t) antitau = sqrt(t*t-antiz*antiz);         // determine tau (z*z) should always be less than (t*t)
        //cout << "OK before the call continue statement..." << endl;
        if (antitau<hydroTau0 || antitau>hydroTauMax-1) continue;       // do not evolve if tau<tau0 or tau>tauMax
        //cout << "x = " << x << ", y = " << y << ", z = " << z << ", t = " << t << endl;
        //cout << "antix = " << antix << ", antiy = " << antiy << ", antiz = " << antiz << ", t = " << t << endl;

        antihydroInfo = hydroSetup->getHydroValues(antix, antiy, antiz, t);// hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
//               hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
        
        if (!trackHistory)
    {
      // warn if out of range and move on to next parton (should not happen!)
      if ( antiix < 0 || antiix >= ixmax ) 
        {
          cout << "WARNING - x out of range" << endl;
          plist[0]->at(i).frozen(1);
          plist[0]->at(antiI).frozen(1);
          continue;
        }
      if ( antiiy < 0 || antiiy >= ixmax )
        {
          cout << "WARNING - y out of range" << endl;
          plist[0]->at(i).frozen(1);
          plist[0]->at(antiI).frozen(1);
          continue;
      }
      if ( (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) && ( antiiz < 0 || antiiz >= izmax ) ) 
        {
          cout << "WARNING - z out of range" << endl;
          plist[0]->at(i).frozen(1);
          plist[0]->at(antiI).frozen(1);
          continue;
        }
    }
        
        antiT = antihydroInfo.T;
        T_ave = 0.5*(T+antiT);
        
        antivx = antihydroInfo.vx;
        antivy = antihydroInfo.vy;
        antivz = antihydroInfo.vz;
        
        antibeta = sqrt(antivx*antivx+antivy*antivy+antivz*antivz);                           // absolute value of flow velocity in units of c
        antigamma = 1./sqrt(1.-antibeta*antibeta);                            // gamma factor
        antiumu[0] = antigamma;
        antiumu[1] = antigamma*antivx;
        antiumu[2] = antigamma*antivy;
        antiumu[3] = antigamma*antivz;
        
        LorentzBooster(antiumu, antipmu, antipmuRest);                      // Boost into the fluid's rest frame
        antidtflow = dt*antipmuRest[0]/antipmu[0];                          // dtfm in the fluid's rest frame
        
      }
      else // the fixed temperature case:
        {
    antibeta = 0.;
    antigamma = 1.; 
    antiT = fixedT; T_ave = fixedT;
    
    antiumu[0] = 1.; antiumu[1] = antiumu[2] = antiumu[3] = 0.;
    LorentzBooster(antiumu, antipmu, antipmuRest);
    antidtflow = dt*antipmuRest[0]/antipmu[0];
    antix = plist[0]->at(antiI).x();                                  // x value of position in [fm]
    antiy = plist[0]->at(antiI).y();                                  // y value of position in [fm]
    antiz = plist[0]->at(antiI).z();                                  // z value of position in [fm]
        }
      // We sample two triplets of kicks now:
      
      //Finally, we take into account the heavy quark-antiquark interaction. This is done in the lab frame, 
      //and is correct strictly non-relativistically. This is good for observables such as TOTAL yields, and 
      //relative yields of excited states, but is probably very bad for something like a differential yield of 
      //J/psi up to high p_T
      if(T_ave > T_C_HQ){
        plist[0]->at(i).tFinal(t); // update the time - finally will be the final time
        plist[0]->at(antiI).tFinal(t); // update the time - finally will be the final time
        
        plist[0]->at(i).x(x+vecp.px()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
        plist[0]->at(i).y(y+vecp.py()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
        plist[0]->at(i).z(z+vecp.pz()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
        plist[0]->at(antiI).x(antix+antivecp.px()/sqrt(antivecp.pAbs()*antivecp.pAbs()+antiM*antiM)*dtfm);
        plist[0]->at(antiI).y(antiy+antivecp.py()/sqrt(antivecp.pAbs()*antivecp.pAbs()+antiM*antiM)*dtfm);
        plist[0]->at(antiI).z(antiz+antivecp.pz()/sqrt(antivecp.pAbs()*antivecp.pAbs()+antiM*antiM)*dtfm);
        
        double TE = T*sqrt(M*M + pmuRest[1]*pmuRest[1] + pmuRest[2]*pmuRest[2] + pmuRest[3]*pmuRest[3]);
        double pHQ = sqrt(pmuRest[1]*pmuRest[1] + pmuRest[2]*pmuRest[2] + pmuRest[3]*pmuRest[3]);
        if(HQ_rate_type == 0)
        {
          kappa_T = 4.*M_PI*T*T*T/TwoPiTD_HQ;
          eta = 2.*M_PI*T*T/(TwoPiTD_HQ*M);
          kappa_L = kappa_T;
        } else if (HQ_rate_type == 1)
        {
          double A0, B0, B1;
          hqRates->getTandPdependent_coeff(T, 0., A0, B0, B1);
          eta = A0;
          kappa_T = 2.*B0;
          kappa_L = kappa_T;
        	//kappa = 2.*eta*T*M;
        } else if (HQ_rate_type == 2)
        {
          double A0, B0, B1;
          hqRates->getTandPdependent_coeff(T, pHQ, A0, B0, B1);
          kappa_T = 2.*B0;
          eta = A0;
          kappa_L = kappa_T;
        } else if (HQ_rate_type == 3)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (10.87937*T - 0.79860)*HQ_Lattice_rate_factor;
          } else
          {
            temp_2piTD = (10.87937*T - 0.79860)*HQ_Lattice_rate_factor;
          }
          kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          //kappa_T = 4.*4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 4)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (10.87937*T - 0.79860);
          } else
          {
            temp_2piTD = 10.87937*T - 0.79860;
          }
          kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
	  eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
	  kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 5)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (6.07614*T - 0.19409);
          } else
          {
            temp_2piTD = 6.07614*T - 0.19409;
          }
        	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
          kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 6)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (15.682606*T - 1.40311);
          } else
          {
            temp_2piTD = 15.682606*T - 1.40311;
          }
        	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
          kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else
        {
        	cout << "Invalid HQ_rate_type" << endl;
          exit(1);
        }

        double kicks[3];
        sigma = sqrt(kappa_T*dtflow);
        kicks[0] = gsl_ran_gaussian(gsl_rand, sigma);
        kicks[1] = gsl_ran_gaussian(gsl_rand, sigma);
        sigma = sqrt(kappa_L*dtflow);
        kicks[2] = gsl_ran_gaussian(gsl_rand, sigma);

        double transv1[3], transv2[3];
        if (pmuRest[3] != 0.)
        {
           transv1[0] = 1.;
           transv1[1] = 1.;
           transv1[2] = -1.*(pmuRest[1]+pmuRest[2])/pmuRest[3];
           double normalization = sqrt(transv1[0]*transv1[0] + transv1[1]*transv1[1] + transv1[2]*transv1[2]);
           for (int ii = 0; ii < 3; ii++) transv1[ii] /= normalization;
        } else
        {
           transv1[0] = 0.;
           transv1[1] = 0.;
           transv1[2] = 1.;
        }

        transv2[0] = pmuRest[2]*transv1[2]/pHQ - pmuRest[3]*transv1[1]/pHQ;
        transv2[1] = pmuRest[3]*transv1[0]/pHQ - pmuRest[1]*transv1[2]/pHQ;
        transv2[2] = pmuRest[1]*transv1[1]/pHQ - pmuRest[2]*transv1[0]/pHQ;

        kicktriplet[0] = kicks[0]*transv1[0] + kicks[1]*transv2[0] + kicks[2]*pmuRest[1]/pHQ;
        kicktriplet[1] = kicks[0]*transv1[1] + kicks[1]*transv2[1] + kicks[2]*pmuRest[2]/pHQ;
        kicktriplet[2] = kicks[0]*transv1[2] + kicks[1]*transv2[2] + kicks[2]*pmuRest[3]/pHQ;

        //Update the momentum in the rest frame:
        pmuRestNew[1] = pmuRest[1]*(1.-eta*dtflow)+kicktriplet[0];
        pmuRestNew[2] = pmuRest[2]*(1.-eta*dtflow)+kicktriplet[1];
        pmuRestNew[3] = pmuRest[3]*(1.-eta*dtflow)+kicktriplet[2];
        pmuRestNew[0] = sqrt(M*M+pmuRestNew[1]*pmuRestNew[1]
           +pmuRestNew[2]*pmuRestNew[2]+pmuRestNew[3]*pmuRestNew[3]);
        
        //The flow velocity for boosting back:
        umub[0] = umu[0];
        umub[1] = -umu[1];
        umub[2] = -umu[2];
        umub[3] = -umu[3];
        
        //The boost back into the lab frame:
        LorentzBooster(umub, pmuRestNew, pmuNew);
        
        //Now, the anti-quark:
        double antiTE = antiT*sqrt(M*M + antipmuRest[1]*antipmuRest[1] + antipmuRest[2]*antipmuRest[2] + antipmuRest[3]*antipmuRest[3]);
        double antipHQ = sqrt(antipmuRest[1]*antipmuRest[1] + antipmuRest[2]*antipmuRest[2] + antipmuRest[3]*antipmuRest[3]);
        if(HQ_rate_type == 0)
        {
          kappa_T = 4.*M_PI*antiT*antiT*antiT/TwoPiTD_HQ;
          eta = 0.5*kappa_T/(M*antiT);
          kappa_L = kappa_T;
        } else if (HQ_rate_type == 1)
        {
          double A0, B0, B1;
          hqRates->getTandPdependent_coeff(antiT, 0., A0, B0, B1);
          kappa_T = 2.*B0;
          eta = A0;
          kappa_L = kappa_T;
      	  //kappa = 2.*eta*T*M;
        } else if (HQ_rate_type == 2)
        {
          double A0, B0, B1;
          hqRates->getTandPdependent_coeff(antiT, antipHQ, A0, B0, B1);
          kappa_T = 2.*B0;
          eta = A0;
          kappa_L = kappa_T;
        } else if (HQ_rate_type == 3)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (10.87937*antiT - 0.79860)*HQ_Lattice_rate_factor;
          } else
          {
            temp_2piTD = (10.87937*antiT - 0.79860)*HQ_Lattice_rate_factor;
          }
          kappa_T = 4.*M_PI*antiT*antiT*antiT/temp_2piTD;
          //kappa_T = 4.*4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*antiT*antiT/(temp_2piTD*M);
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 4)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (10.87937*antiT - 0.79860);
          } else
          {
            temp_2piTD = 10.87937*antiT - 0.79860;
          }
          kappa_T = 4.*M_PI*antiT*antiT*antiT/temp_2piTD;
          eta = 2.*M_PI*antiT*antiT/(temp_2piTD*M);
	  eta     *= 1./(1.+0.17*pow(antipHQ,1.11));//pHQ in GeV
	  kappa_T *= (1.+ 5.*pow(antipHQ,1.11)/(7.9+pow(antipHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 5)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (6.07614*T - 0.19409);
          } else
          {
            temp_2piTD = 6.07614*T - 0.19409;
          }
        	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
          kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 6)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (15.682606*T - 1.40311);
          } else
          {
            temp_2piTD = 15.682606*T - 1.40311;
          }
        	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
          kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else
        {
        	cout << "Invalid HQ_rate_type" << endl;
          exit(1);
        }
        sigma = sqrt(kappa_T*dtflow);
        kicks[0] = gsl_ran_gaussian(gsl_rand, sigma);
        kicks[1] = gsl_ran_gaussian(gsl_rand, sigma);
        sigma = sqrt(kappa_L*dtflow);
        kicks[2] = gsl_ran_gaussian(gsl_rand, sigma);

        double antitransv1[3], antitransv2[3];
        if (antipmuRest[3] != 0.)
        {
           antitransv1[0] = 1.;
           antitransv1[1] = 1.;
           antitransv1[2] = -1.*(antipmuRest[1]+antipmuRest[2])/antipmuRest[3];
           double normalization = sqrt(antitransv1[0]*antitransv1[0] + antitransv1[1]*antitransv1[1] + antitransv1[2]*antitransv1[2]);
           for (int ii = 0; ii < 3; ii++) antitransv1[ii] /= normalization;
        } else
        {
           antitransv1[0] = 0.;
           antitransv1[1] = 0.;
           antitransv1[2] = 1.;
        }

        antitransv2[0] = antipmuRest[2]*antitransv1[2]/antipHQ - antipmuRest[3]*antitransv1[1]/antipHQ;
        antitransv2[1] = antipmuRest[3]*antitransv1[0]/antipHQ - antipmuRest[1]*antitransv1[2]/antipHQ;
        antitransv2[2] = antipmuRest[1]*antitransv1[1]/antipHQ - antipmuRest[2]*antitransv1[0]/antipHQ;

        kicktriplet[0] = kicks[0]*antitransv1[0] + kicks[1]*antitransv2[0] + kicks[2]*antipmuRest[1]/antipHQ;
        kicktriplet[1] = kicks[0]*antitransv1[1] + kicks[1]*antitransv2[1] + kicks[2]*antipmuRest[2]/antipHQ;
        kicktriplet[2] = kicks[0]*antitransv1[2] + kicks[1]*antitransv2[2] + kicks[2]*antipmuRest[3]/antipHQ;
        
        //Update the momentum in the rest frame:
        antipmuRestNew[1] = antipmuRest[1]*(1.-eta*antidtflow)+kicktriplet[0];
        antipmuRestNew[2] = antipmuRest[2]*(1.-eta*antidtflow)+kicktriplet[1];
        antipmuRestNew[3] = antipmuRest[3]*(1.-eta*antidtflow)+kicktriplet[2];
        antipmuRestNew[0] = sqrt(M*M+antipmuRestNew[1]*antipmuRestNew[1]
               +antipmuRestNew[2]*antipmuRestNew[2]+antipmuRestNew[3]*antipmuRestNew[3]);
        
        //The flow velocity for boosting back:
        antiumub[0] = antiumu[0];
        antiumub[1] = -antiumu[1];
        antiumub[2] = -antiumu[2];
        antiumub[3] = -antiumu[3];
        
        //The boost back into the lab frame:
        LorentzBooster(antiumub, antipmuRestNew, antipmuNew);
        
        //Finally, we include the effect of the HQ potential, parametrized from
        // lattice data:
        r = sqrt( (x-antix)*(x-antix)+(y-antiy)*(y-antiy)+(z-antiz)*(z-antiz) );
        force = F(r, T_ave);
        //force = 0.; //We are debugging!
        //if(r > 0.4 && r < 1.5){
        //force = 0.0985/(r*r)+0.812;
        //}//Still debugging!

        if(force > 0.){
    F_x = force*(antix-x)/r;
    F_y = force*(antiy-y)/r;
    F_z = force*(antiz-z)/r;
        }
        else F_x = F_y = F_z = 0.;

        //Update the particle's momentum:
        plist[0]->at(i).p(pmuNew[1]+F_x*dt*hbarc, pmuNew[2]+F_y*dt*hbarc, pmuNew[3]+F_z*dt*hbarc);
        plist[0]->at(antiI).p(antipmuNew[1]-F_x*dt*hbarc, antipmuNew[2]-F_y*dt*hbarc, antipmuNew[3]-F_z*dt*hbarc);
        
        //Update the anti-particle's momentum:
        //cout << "pmu1: " << pmuNew[1] << " " << pmuNew[2] << " " << pmuNew[3] << endl;
        //cout << "pmu2: " << antipmuNew[1] << " " << antipmuNew[2] << " " << antipmuNew[3] << endl;
      }
      else{
        plist[0]->at(i).frozen(1);
        plist[0]->at(antiI).frozen(1);
      }
    }
  }
      }
    }
  //foutt.close();
  return counter;
}

//Evolves the heavy quarks, but does the hadronization each heavy quark reaches the temperature of 
//its freezeout. This is to make better estimates of the recombinant production. -CFY 8/23/2011
int MARTINI::evolveAndHadronizeHeavyQuarks(vector<Parton>  **plist, double* hydro_T, double* prehydro_T, int* hydro_counter, int* prehydro_counter, int counter, int it)
{
  cout.precision(6);
  
  HydroInfo hydroInfo;// antihydroInfo;
  PreHydroInfo prehydroInfo;// antihydroInfo;
  
  Vec4 vecp, antiVecP;                          // momentum four-vectors
  double pmu[4], pmuRest[4], pmuRestNew[4], pmuNew[4]; //Easy to Lorentz-boost arrays with my subroutine 
  //double antipmu[4], antipmuRest[4], antipmuRestNew[4], antipmuNew[4]; //Easy to Lorentz-boost arrays with my subroutine 
  
  int imax = plist[0]->size();        // number of partons in the list
  int id; //newID;                      // original and new parton's ID
  //int col, acol, newCol, newAcol;     // number of the color string attached to the old, and new parton
  int ix, iy, iz, itau;               // parton's position in cell coordinates
  //int antiix, antiiy, antiiz, antiitau; // parton's position in cell coordinates
  int ixmax, izmax;                   // maximum cell number in x- and y-direction, and z-direction
  
  double dt = dtfm/hbarc;             // time step in [GeV^(-1)] - hbarc is defined in Basics.h 
  // Gamma = \int dGamma dp is given in units of GeV, so dt*Gamma is dimensionless
  double dtflow;// antidtflow;                      // dt in the fluid's rest frame
  double T;// antiT, T_ave;             // temperature
  double M, antiM;                    // Parton's mass [GeV]
  double x, y, z;                     // parton's position in [fm]
  //double antix, antiy, antiz       ;  // parton's position in [fm]
  double t, tau;// antitau;             // lab time and lab tau
  double vx, vy, vz;                  // flow velocity of the cell
  //double antivx, antivy, antivz;      // flow velocity of the cell
  double beta;                        // absolute value of the flow velocity
  //double antibeta;                        // absolute value of the flow velocity
  double gamma;                       // gamma factor 
  //double antigamma;                       // gamma factor 
  double umu[4], umub[4];             // Fluid's 4-velocity
  //double antiumu[4], antiumub[4];             // Fluid's 4-velocity
  const double rd=12.3;               // ratio of degrees of freedom g_QGP/g_H (for Kolb hydro only)
  
  double kappa_T, kappa_L, sigma, eta;
  double kicktriplet[3];
  double r, force, F_x, F_y, F_z;
  
  ixmax = static_cast<int>(2.*hydroXmax/hydroDx+0.0001);
  
  if ( hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6) izmax = static_cast<int>(2.*hydroZmax/hydroDz+0.0001);
  
  //cout << "counter = " << counter << ", it = " << it << endl;
  
  for ( int i=0; i<imax; i++)                                   // loop over all partons
    { 
      if(plist[0]->at(i).frozen() == 0)
      {
        counter++;                                                // counter that characterizes every step uniquely
        id = plist[0]->at(i).id();                                // id of parton i 
        
        if(abs(id)!=4 && abs(id)!=5) continue;                    // We immediately leave the loop of the parton is not a heavy quark!
        
        
        t = it*dtfm;        //had a +tau0 here. now let partons travel doing their vacuum shower stuff first...
        
        vecp = plist[0]->at(i).p();                               // four-vector p of parton i
        M = plist[0]->at(i).mass();                               // Mass of the parton (only used for charm and bottom quarks)
        pmu[1] = vecp.px();
        pmu[2] = vecp.py();
        pmu[3] = vecp.pz();
        pmu[0] = sqrt(M*M+pmu[1]*pmu[1]+pmu[2]*pmu[2]+pmu[3]*pmu[3]); //vecp loaded into an array for ease with LorentzBooster
        
	if(plist[0]->at(i).formed() == 0)
        {
          if (pmu[0]*hbarc/(M*M) < t) plist[0]->at(i).formed(1);
        }

        x = plist[0]->at(i).x();                                  // x value of position in [fm]
        y = plist[0]->at(i).y();                                  // y value of position in [fm]
        z = plist[0]->at(i).z();                                  // z value of position in [fm]
        
        ix = floor((hydroXmax+x)/hydroDx+0.0001);                 // x-coordinate of the cell we are in now
        iy = floor((hydroXmax+y)/hydroDx+0.0001);                 // y-coordinate of the cell we are in now
        
        if (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) 
        {
          iz = floor((hydroZmax+z)/hydroDz+0.0001);             // z-coordinate of the cell we are in now
        }
        
        if (fixedTemperature == 0)
        {
          tau = 0.;           //just give tau some value - will change below
          
          if (z*z<t*t) tau = sqrt(t*t-z*z);                     // determine tau (z*z) should always be less than (t*t)
	  if(tau >= hydroTau0 && plist[0]->at(i).inhydro() == 0)
          {
            Vec4 initp;
	    initp = plist[0]->at(i).Initialp();
            plist[0]->at(i).EndOfPreeqp(pmu[1], pmu[2], pmu[3]);
            plist[0]->at(i).inhydro(1);
          }
          // get the temperature and flow velocity at the current position and time:
          if (!trackHistory)
            {
              // warn if out of range and move on to next parton (should not happen!)
              if ( ix < 0 || ix >= ixmax ) 
              {
                cout << "WARNING - x out of range" << endl;
                plist[0]->at(i).frozen(1);
                continue;
              }
              if ( iy < 0 || iy >= ixmax )
              {
                cout << "WARNING - y out of range" << endl;
                plist[0]->at(i).frozen(1);
                continue;
              }
              if ( (hydroWhichHydro == 3 || hydroWhichHydro == 4 || hydroWhichHydro == 6 ) && ( iz < 0 || iz >= izmax ) ) 
              {
                cout << "WARNING - z out of range" << endl;
                plist[0]->at(i).frozen(1);
                continue;
              }
            }
          if(evolution == 1)
	    {
              if(tau<hydroTau0 && tau > prehydroTau0 && prehydroevolution == 1 && plist[0]->at(i).formed() == 1)
              {
                prehydroInfo = hydroSetup->getPreHydroValues(x, y, z, t);//, hydroXmax, hydroZmax, hydroTauMax, hydroTau0
                T = prehydroInfo.T;
                double eta_local = 0.5*log((t+z)/(t-z));
		double cosh_eta = cosh(eta_local);
		double sinh_eta = sinh(eta_local);
		umu[0] = prehydroInfo.utau*cosh_eta;  // gamma factor
                umu[1] = prehydroInfo.ux; 
                umu[2] = prehydroInfo.uy; 
                umu[3] = prehydroInfo.utau*sinh_eta;

		LorentzBooster(umu, pmu, pmuRest);                      // Boost into the fluid's rest frame
                dtflow = dt*pmuRest[0]/pmu[0];                          // dtfm in the fluid's rest frame
                if (T > T_C_HQ){
                  int itind = (tau - prehydroTau0)/prehydroDtau;
                  prehydro_T[itind] += T;
                  prehydro_counter[itind] ++;
                }
              } else if (tau >= hydroTau0 && plist[0]->at(i).formed() == 1)
              {
                hydroInfo = hydroSetup->getHydroValues(x, y, z, t);//, hydroXmax, hydroZmax, hydroTauMax, hydroTau0, 
//                       hydroDx, hydroDz, hydroDtau, hydroWhichHydro, fixedDistribution, lattice, trackHistory);
              
                T = hydroInfo.T;
                //cout << "T = " << T << "for particle " << i << endl;
                vx = hydroInfo.vx;
                vy = hydroInfo.vy;
                vz = hydroInfo.vz;
                beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of flow velocity in units of c
                gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
                umu[0] = gamma;
                umu[1] = gamma*vx;
                umu[2] = gamma*vy;
                umu[3] = gamma*vz;
                LorentzBooster(umu, pmu, pmuRest);                      // Boost into the fluid's rest frame
                dtflow = dt*pmuRest[0]/pmu[0];                          // dtfm in the fluid's rest frame
                if (T > T_C_HQ){
                  int itind = (tau - hydroTau0)/hydroDtau;
                  hydro_T[itind] += T;
                  hydro_counter[itind] ++;
                }
                
              }
              else
              {
                T = 0.;
              }
            }
          else{
              //This is how to stop evolution in medium if the settings call for it:
              T = 0.;
            }

	  //We add here a third option for moveBeforeTau0: if moveBeforeTau0 = 2, then the heavy quark will be evolved according to the 
          //Cornell potential if tau<hydroTau0:
          if(plist[0]->at(i).formed() == 0 || (tau<hydroTau0 && moveBeforeTau0 == 2 && (prehydroevolution == 0 || T < T_C_HQ)))
          {
            //First, evolve the positions:
            plist[0]->at(i).x(x+vecp.px()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
            plist[0]->at(i).y(y+vecp.py()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
            plist[0]->at(i).z(z+vecp.pz()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);

            plist[0]->at(i).xOld(x);
            plist[0]->at(i).yOld(y);
            plist[0]->at(i).zOld(z);

	    //MS Turning off Cornell potential
            //int antiI = plist[0]->at(i).antiI();
            ////cout << "antiI = " << antiI << " for i = " << i << endl;
            //if(plist[0]->at(antiI).frozen() == 0)
            //{
            //  //Here are the steps necessary for evolving in the CM frame of the pair:
            //  if(antiI > i){
            //   antiVecP = plist[0]->at(antiI).p();
            //  }
            //  else{
            //   antiVecP = plist[0]->at(antiI).pOld();
            //  }
            //  antiM = plist[0]->at(antiI).mass();
            //  
            //  double dx, dy, dz;
            //  if(antiI > i)
            //  {
            //    dx = plist[0]->at(antiI).x()+antiVecP.px()/sqrt(antiVecP.pAbs()*antiVecP.pAbs()+antiM*antiM)*dtfm-plist[0]->at(i).x();
            //    dy = plist[0]->at(antiI).y()+antiVecP.py()/sqrt(antiVecP.pAbs()*antiVecP.pAbs()+antiM*antiM)*dtfm-plist[0]->at(i).y();
            //    dz = plist[0]->at(antiI).z()+antiVecP.pz()/sqrt(antiVecP.pAbs()*antiVecP.pAbs()+antiM*antiM)*dtfm-plist[0]->at(i).z();
            //  }
            //  else{
            //    dx = plist[0]->at(antiI).x()-plist[0]->at(i).x();
            //    dy = plist[0]->at(antiI).y()-plist[0]->at(i).y();
            //    dz = plist[0]->at(antiI).z()-plist[0]->at(i).z();
            //  }
            //  
            //  int id1 = plist[0]->at(i).id();
            //  int id2 = plist[0]->at(antiI).id();
            //  
            //  double ptot[4];
            //  ptot[1] = vecp.px() + antiVecP.px();
            //  ptot[2] = vecp.py() + antiVecP.py();
            //  ptot[3] = vecp.pz() + antiVecP.pz();
            //  ptot[0] = sqrt(vecp.px()*vecp.px() + vecp.py()*vecp.py() + vecp.pz()*vecp.pz() + M*M)
            //     +sqrt(antiVecP.px()*antiVecP.px() + antiVecP.py()*antiVecP.py() + antiVecP.pz()*antiVecP.pz() + antiM*antiM);
            //  
            //  //The space-like separation between the heavy quarks:
            //  double deltax[4];
            //  deltax[0] = 0.;
            //  deltax[1] = dx;
            //  deltax[2] = dy;
            //  deltax[3] = dz;
            //  
            //  double pnew[4];
            //  
            //  //Now we have everything needed to calculate pnew:
            //  stepMomentumInCMFrame(ptot, pmu, dt, deltax, pnew, 0.);
            //  
            //  double px = vecp.px();
            //  double py = vecp.py();
            //  double pz = vecp.pz();
            //  
            //  plist[0]->at(i).p(pnew[1], pnew[2], pnew[3]);
            //  plist[0]->at(i).pOld(px, py, pz);
            //}
          }

          if ((tau<hydroTau0 && prehydroevolution == 0) || (tau>hydroTauMax && evolution == 1) || tau < prehydroTau0 || plist[0]->at(i).formed() == 0) continue;       // do not evolve if tau<tau0 or tau>tauMax
          
          //This is how to stop evolution in medium if the settings call for it:


          if (T<hydroTfinal && hydroWhichHydro == 2)                 // for 2D hydro evolutions stop at T_final (can change that)
            {
              plist[0]->at(i).tFinal(t);
              continue;                                             // do not evolve if T<T_c
            }
          
        }
            else // the fixed temperature case:
        {
          beta = 0.;
          gamma = 1.; 
          T = fixedT;
          umu[0] = 1.; umu[1] = umu[2] = umu[3] = 0.;
          LorentzBooster(umu, pmu, pmuRest);
          dtflow = dt*pmuRest[0]/pmu[0];
        }
        // Sample a triplet of kicks for the heavy quark. For now, we 
        // take the simplest case: constant 2piT*D_c, however later 
        // both kicktriplet[3] and eta will be determined from HTL 
        // results at LO or NLO order:
    if(T > T_C_HQ){
      plist[0]->at(i).tFinal(t); // update the time - finally will be the final time

      //Finally evolve the heavy quarks' positions:
      plist[0]->at(i).x(x+vecp.px()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
      plist[0]->at(i).y(y+vecp.py()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);
      plist[0]->at(i).z(z+vecp.pz()/sqrt(vecp.pAbs()*vecp.pAbs()+M*M)*dtfm);

      plist[0]->at(i).xOld(x);
      plist[0]->at(i).yOld(y);
      plist[0]->at(i).zOld(z);
      
      int antiI = plist[0]->at(i).antiI();

      if(plist[0]->at(antiI).frozen() == 0){
        //Here are the steps necessary for evolving in the CM frame of the pair:
        if(antiI > i){
          antiVecP = plist[0]->at(antiI).p();
        }
        else{
          antiVecP = plist[0]->at(antiI).pOld();
        }
        antiM = plist[0]->at(antiI).mass();
        
        double dx, dy, dz;
        if(antiI > i){
          dx = plist[0]->at(antiI).x()+antiVecP.px()/sqrt(antiVecP.pAbs()*antiVecP.pAbs()+antiM*antiM)*dtfm-plist[0]->at(i).x();
          dy = plist[0]->at(antiI).y()+antiVecP.py()/sqrt(antiVecP.pAbs()*antiVecP.pAbs()+antiM*antiM)*dtfm-plist[0]->at(i).y();
          dz = plist[0]->at(antiI).z()+antiVecP.pz()/sqrt(antiVecP.pAbs()*antiVecP.pAbs()+antiM*antiM)*dtfm-plist[0]->at(i).z();
        }
        else{
          dx = plist[0]->at(antiI).x()-plist[0]->at(i).x();
          dy = plist[0]->at(antiI).y()-plist[0]->at(i).y();
          dz = plist[0]->at(antiI).z()-plist[0]->at(i).z();
        }
        
        int id1 = plist[0]->at(i).id();
        int id2 = plist[0]->at(antiI).id();
        
        double ptot[4];
        ptot[1] = vecp.px() + antiVecP.px();
        ptot[2] = vecp.py() + antiVecP.py();
        ptot[3] = vecp.pz() + antiVecP.pz();
        ptot[0] = sqrt(vecp.px()*vecp.px() + vecp.py()*vecp.py() + vecp.pz()*vecp.pz() + M*M)
          +sqrt(antiVecP.px()*antiVecP.px() + antiVecP.py()*antiVecP.py() + antiVecP.pz()*antiVecP.pz() + antiM*antiM);
        
        //The space-like separation between the heavy quarks:
        double deltax[4];
        deltax[0] = 0.;
        deltax[1] = dx;
        deltax[2] = dy;
        deltax[3] = dz;
        
        double pnew[4];
        
        //Now we have everything needed to calculate pnew:
        stepMomentumInCMFrame(ptot, pmu, dt, deltax, pnew, T);
        
        double px = vecp.px();
        double py = vecp.py();
        double pz = vecp.pz();
        
        plist[0]->at(i).p(pnew[1], pnew[2], pnew[3]);
        plist[0]->at(i).pOld(px, py, pz);
        vecp = plist[0]->at(i).p();
      }

      //Reload and recalculate these for the Langevin dynamics:
      pmu[1] = vecp.px();
      pmu[2] = vecp.py();
      pmu[3] = vecp.pz();
      pmu[0] = sqrt(M*M+pmu[1]*pmu[1]+pmu[2]*pmu[2]+pmu[3]*pmu[3]); //vecp loaded into an array for ease with LorentzBooster
      LorentzBooster(umu, pmu, pmuRest);                      // Boost into the fluid's rest frame
      dtflow = dt*pmuRest[0]/pmu[0];                          // dtfm in the fluid's rest frame

      double TE = T*sqrt(M*M + pmuRest[1]*pmuRest[1] + pmuRest[2]*pmuRest[2] + pmuRest[3]*pmuRest[3]);
      double pHQ = sqrt(pmuRest[1]*pmuRest[1] + pmuRest[2]*pmuRest[2] + pmuRest[3]*pmuRest[3]);
      if(HQ_rate_type == 0)
      {
        kappa_T = 4.*M_PI*T*T*T/TwoPiTD_HQ;
        eta = 2.*M_PI*T*T/(TwoPiTD_HQ*M);
        kappa_L = kappa_T;
      } else if (HQ_rate_type == 1)
      {
        double A0, B0, B1;
        hqRates->getTandPdependent_coeff(T, 0., A0, B0, B1);
        eta = A0;
        kappa_T = 2.*B0;
        kappa_L = kappa_T;
      	//kappa = 2.*eta*T*M;
      } else if (HQ_rate_type == 2)
      {
        double A0, B0, B1;
        hqRates->getTandPdependent_coeff(T, pHQ, A0, B0, B1);
      	kappa_T = 2.*B0;
        eta = A0;
        kappa_L = kappa_T;
      } else if (HQ_rate_type == 3)
      {
        double temp_2piTD;
        if (tau < hydroTau0)
        {
//          temp_2piTD = 20.30405*T + 0.35256;
          temp_2piTD = (10.87937*T - 0.79860)*HQ_Lattice_rate_factor;
        } else
        {
          temp_2piTD = (10.87937*T - 0.79860)*HQ_Lattice_rate_factor;
        }
        kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
        //kappa_T = 4.*4.*M_PI*T*T*T/temp_2piTD;
        eta = 2.*M_PI*T*T/(temp_2piTD*M);
        kappa_L = kappa_T;

        } else if (HQ_rate_type == 4)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (10.87937*T - 0.79860);
          } else
          {
            temp_2piTD = 10.87937*T - 0.79860;
          }
        	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
          kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 5)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (6.07614*T - 0.19409);
          } else
          {
            temp_2piTD = 6.07614*T - 0.19409;
          }
        	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
          kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 6)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (15.682606*T - 1.40311);
          } else
          {
            temp_2piTD = 15.682606*T - 1.40311;
          }
        	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          eta     *= 1./(1.+0.17*pow(pHQ,1.11));//pHQ in GeV
          kappa_T *= (1.+ 5.*pow(pHQ,1.11)/(7.9+pow(pHQ,1.11)));//pHQ in GeV
          kappa_L = kappa_T;

        } else if (HQ_rate_type == 7)
        {
          double temp_2piTD;
          if (tau < hydroTau0)
          {
//            temp_2piTD = 20.30405*T + 0.35256;
            temp_2piTD = (10.87937*T - 0.79860);
          } else
          {
            temp_2piTD = 10.87937*T - 0.79860;
          }
        	kappa_T = 4.*M_PI*T*T*T/temp_2piTD;
          eta = 2.*M_PI*T*T/(temp_2piTD*M);
          //eta     *= 1./(1.+0.1*pHQ);//pHQ in GeV
          eta     *= 1./(1.+0.2*pow(pHQ,0.5));//pHQ in GeV
          kappa_T *= (1.+ 5.*pHQ/(8+pHQ));//pHQ in GeV
          kappa_L = kappa_T;

        } else
        {
        	cout << "Invalid HQ_rate_type" << endl;
          exit(1);
        }

      double kicks[3];
      sigma = sqrt(kappa_T*dtflow);
      kicks[0] = gsl_ran_gaussian(gsl_rand, sigma);
      kicks[1] = gsl_ran_gaussian(gsl_rand, sigma);
      sigma = sqrt(kappa_L*dtflow);
      kicks[2] = gsl_ran_gaussian(gsl_rand, sigma);

      double transv1[3], transv2[3];
      if (pmuRest[3] != 0.)
      {
         transv1[0] = 1.;
         transv1[1] = 1.;
         transv1[2] = -1.*(pmuRest[1]+pmuRest[2])/pmuRest[3];
         double normalization = sqrt(transv1[0]*transv1[0] + transv1[1]*transv1[1] + transv1[2]*transv1[2]);
         for (int ii = 0; ii < 3; ii++) transv1[ii] /= normalization;
      } else
      {
         transv1[0] = 0.;
         transv1[1] = 0.;
         transv1[2] = 1.;
      }

      transv2[0] = pmuRest[2]*transv1[2]/pHQ - pmuRest[3]*transv1[1]/pHQ;
      transv2[1] = pmuRest[3]*transv1[0]/pHQ - pmuRest[1]*transv1[2]/pHQ;
      transv2[2] = pmuRest[1]*transv1[1]/pHQ - pmuRest[2]*transv1[0]/pHQ;

      kicktriplet[0] = kicks[0]*transv1[0] + kicks[1]*transv2[0] + kicks[2]*pmuRest[1]/pHQ;
      kicktriplet[1] = kicks[0]*transv1[1] + kicks[1]*transv2[1] + kicks[2]*pmuRest[2]/pHQ;
      kicktriplet[2] = kicks[0]*transv1[2] + kicks[1]*transv2[2] + kicks[2]*pmuRest[3]/pHQ;
      
      //Update the momentum in the rest frame:
      pmuRestNew[1] = pmuRest[1]*(1.-eta*dtflow)+kicktriplet[0];
      pmuRestNew[2] = pmuRest[2]*(1.-eta*dtflow)+kicktriplet[1];
      pmuRestNew[3] = pmuRest[3]*(1.-eta*dtflow)+kicktriplet[2];
      pmuRestNew[0] = sqrt(M*M+pmuRestNew[1]*pmuRestNew[1]
               +pmuRestNew[2]*pmuRestNew[2]+pmuRestNew[3]*pmuRestNew[3]);

//      if (tau < hydroTau0) {      
//        pmuRestNew[1] = pmuRest[1];//*(1.-eta*dtflow)+kicktriplet[0];
//        pmuRestNew[2] = pmuRest[2];//*(1.-eta*dtflow)+kicktriplet[1];
//        pmuRestNew[3] = pmuRest[3]*(1.-eta*dtflow)+kicktriplet[2];
//        pmuRestNew[0] = sqrt(M*M+pmuRestNew[1]*pmuRestNew[1]
//                 +pmuRestNew[2]*pmuRestNew[2]+pmuRestNew[3]*pmuRestNew[3]);
//      } else {
//        pmuRestNew[1] = pmuRest[1]*(1.-eta*dtflow)+kicktriplet[0];
//        pmuRestNew[2] = pmuRest[2]*(1.-eta*dtflow)+kicktriplet[1];
//        pmuRestNew[3] = pmuRest[3]*(1.-eta*dtflow)+kicktriplet[2];
//        pmuRestNew[0] = sqrt(M*M+pmuRestNew[1]*pmuRestNew[1]
//                 +pmuRestNew[2]*pmuRestNew[2]+pmuRestNew[3]*pmuRestNew[3]);
//      }
      //The flow velocity for boosting back:
      umub[0] = umu[0];
      umub[1] = -umu[1];
      umub[2] = -umu[2];
      umub[3] = -umu[3];
      
      //The boost back into the lab frame:
      LorentzBooster(umub, pmuRestNew, pmuNew);
      
      plist[0]->at(i).p(pmuNew[1], pmuNew[2], pmuNew[3]);
      
          }
          // If the heavy quark's temperature is below that of freezeout, hadronize the heavy quark, 
          // taking into account the possibility of recombinant production. Recombinant production is assumed to be 
          // possible only if both heavy quarks involved had evolved in the medium (in other words, if tFinal != 0 fm/c):
          else if (T < T_C_HQ && (prehydroevolution == 0 || tau > hydroTau0) && plist[0]->at(i).formed() == 1)
       {
      //cout << "Is this being called?" << endl;
      plist[0]->at(i).frozen(1);
      //cout << "The heavy quark " << i << " is hadronizing..." << endl;

      //First, check if quarkonium has survived:
      int antiI = plist[0]->at(i).antiI();

      if(plist[0]->at(antiI).frozen() == 0){
        int id1 = plist[0]->at(i).id();
        int id2 = plist[0]->at(antiI).id();

        //Any quarkonium produced will be the surviving component:
        int status = 81;

        //The relative momentum of the pair:
        Vec4 p1 = plist[0]->at(i).p();
        Vec4 p2 = plist[0]->at(antiI).p();

        double M1 = plist[0]->at(i).mass();
        double M2 = plist[0]->at(antiI).mass();

        double p1mu[4], p2mu[4];
        p1mu[1] = p1.px();
        p1mu[2] = p1.py();
        p1mu[3] = p1.pz();
        p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
        p2mu[1] = p2.px();
        p2mu[2] = p2.py();
        p2mu[3] = p2.pz();
        p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);
        
        double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
                -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
                -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
                -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );
        
        double deltaX[4];
        deltaX[0] = 0.;
        deltaX[1] = plist[0]->at(i).x()-plist[0]->at(antiI).x();
        deltaX[2] = plist[0]->at(i).y()-plist[0]->at(antiI).y();
        deltaX[3] = plist[0]->at(i).z()-plist[0]->at(antiI).z();

        double r = rInCMFrame(p1mu, p2mu, deltaX);

        double potential = CornellPotential(r);

        //Finally, the energy of the pair in the CM frame:
        double E_CM = M_inv+potential;

        int bound = isBound(id1, id2, E_CM);

        if(bound == 1){
          //cout << "Surviving quarkonium should be formed, with status = " << status << endl;
          //cout << "tFinal of parton " << i << " = " << plist[0]->at(i).tFinal() << endl;
          //cout << "tFinal of parton " << antiI << " = " << plist[0]->at(antiI).tFinal() << endl;
          hadronizeByColorEvaporation(id1, id2, E_CM, p1mu, p2mu, M1, M2, status);  
          plist[0]->at(antiI).frozen(1);
          //If this led to hadronization, break the loop here:
          continue;
        }
      }

      //If the heavy quark did not evolve in the medium, fragment its momentum and hadronize it to 
      //an open heavy flavor meson:
      if((plist[0]->at(i).tFinal() < hydroTau0 && prehydroevolution == 0) || plist[0]->at(i).tFinal() < prehydroTau0){
        //cout << "Unevolved open heavy flavor should be formed:" << endl;
        //cout << "tFinal of parton " << i << " = " << plist[0]->at(i).tFinal() << endl;
        Vec4 p1 = plist[0]->at(i).p();
        double M1 = plist[0]->at(i).mass();
        int id1 = plist[0]->at(i).id();
        double p1mu[4];
        p1mu[1] = p1.px();
        p1mu[2] = p1.py();
        p1mu[3] = p1.pz();
        p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
        hadronizeOneHeavyQuark(id1, p1mu, M1);
        continue;
      }

      //If the quark did not hadronize by surviving as initially-produced quarkonium, check 
      //for recombinant production. For now, all pairings with neighbors are randomly selected:
      //cout << "Now we are checking for recombinant production..." << endl;
      int id1 = plist[0]->at(i).id();
      int charge1 = id1/abs(id1);
      vector<int> neighbors;
      for(int j=0; j<imax; j++){
        int id2 = plist[0]->at(j).id();
        int charge2 = id2/abs(id2);
        double evoltime;
        if (prehydroevolution == 0)
        {
           evoltime = hydroTau0;
        } else {
           evoltime = prehydroTau0;
        }
        if(charge2 != charge1
           && plist[0]->at(j).tFinal() > evoltime){
          //The relative momentum of the pair:
          Vec4 p1 = plist[0]->at(i).p();
          Vec4 p2 = plist[0]->at(j).p();
          
          double M1 = plist[0]->at(i).mass();
          double M2 = plist[0]->at(j).mass();
          
          double p1mu[4], p2mu[4];
          p1mu[1] = p1.px();
          p1mu[2] = p1.py();
          p1mu[3] = p1.pz();
          p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
          p2mu[1] = p2.px();
          p2mu[2] = p2.py();
          p2mu[3] = p2.pz();
          p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);
          
          double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
            -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
            -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
            -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );
          
          double deltaX[4];
          deltaX[0] = 0.;
          deltaX[1] = plist[0]->at(i).x()-plist[0]->at(j).x();
          deltaX[2] = plist[0]->at(i).y()-plist[0]->at(j).y();
          deltaX[3] = plist[0]->at(i).z()-plist[0]->at(j).z();

          double r = rInCMFrame(p1mu, p2mu, deltaX);
          double potential = CornellPotential(r);
          
          //Finally, the energy of the pair in the CM frame:
          double E_CM = M_inv+potential;
          
          int bound = isBound(id1, id2, E_CM);
          if(bound == 1){
            neighbors.push_back(j);
          }
        }
      }
      //Randomly select a quark that fits the criteria:
      int numberOfNeighbors = neighbors.size();

      if(numberOfNeighbors == 0){
        Vec4 p1 = plist[0]->at(i).p();
        double M1 = plist[0]->at(i).mass();
        double x = plist[0]->at(i).x();                                  // x value of position in [fm]
        double y = plist[0]->at(i).y();                                  // y value of position in [fm]
        double z = plist[0]->at(i).z();                                  // z value of position in [fm]
        double t = plist[0]->at(i).tFinal();
        double p1mu[4];
        p1mu[1] = p1.px();
        p1mu[2] = p1.py();
        p1mu[3] = p1.pz();
        p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
        if (do_coalesence)
        {
            int does_it_coalesce = check_coalescence(id1, p1mu, M1, x, y, z, t);
            if (does_it_coalesce == 1 && id1 == 4) {
                bool complete = hadronize_by_coalescence(id1, p1mu, 0.3, x, y, z, t);
                if (!complete) hadronizeOneHeavyQuark(id1, p1mu, M1);
            } else if (does_it_coalesce == 0 || id1 != 4) {
                hadronizeOneHeavyQuark(id1, p1mu, M1);
            }
        } else {
            hadronizeOneHeavyQuark(id1, p1mu, M1);
        }
        continue;
      }

      int iPicked = gsl_rng_uniform_int(gsl_rand, numberOfNeighbors);
      //cout << "iPicked = " << iPicked << endl;
      if(iPicked == numberOfNeighbors){
        Vec4 p1 = plist[0]->at(i).p();
        double M1 = plist[0]->at(i).mass();
        double p1mu[4];
        p1mu[1] = p1.px();
        p1mu[2] = p1.py();
        p1mu[3] = p1.pz();
        p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
        if (do_coalesence)
        {
            int does_it_coalesce = check_coalescence(id1, p1mu, M1, x, y, z, t);
            if (does_it_coalesce == 1 && id1 == 4 && do_coalesence) {
                bool complete = hadronize_by_coalescence(id1, p1mu, 0.3, x, y, z, t);
                if (!complete) hadronizeOneHeavyQuark(id1, p1mu, M1);
            } else if (does_it_coalesce == 0 || id1 != 4) {
                hadronizeOneHeavyQuark(id1, p1mu, M1);
            }
        } else {
            hadronizeOneHeavyQuark(id1, p1mu, M1);
        }
      }
      else{
        int i2 = neighbors.at(iPicked);
        int id2 = plist[0]->at(i2).id();

        //The relative momentum of the pair:
        Vec4 p1 = plist[0]->at(i).p();
        Vec4 p2 = plist[0]->at(i2).p();
          
        double M1 = plist[0]->at(i).mass();
        double M2 = plist[0]->at(i2).mass();
        
        double p1mu[4], p2mu[4];
        p1mu[1] = p1.px();
        p1mu[2] = p1.py();
        p1mu[3] = p1.pz();
        p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
        p2mu[1] = p2.px();
        p2mu[2] = p2.py();
        p2mu[3] = p2.pz();
        p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);
        
        double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
                -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
                -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
                -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );

        //double dx, dy, dz;
        //dx = plist[0]->at(i).x()-plist[0]->at(i2).x();
        //dy = plist[0]->at(i).y()-plist[0]->at(i2).y();
        //dz = plist[0]->at(i).z()-plist[0]->at(i2).z();
        //double r = sqrt(dx*dx+dy*dy+dz*dz);

        double deltaX[4];
        deltaX[0] = 0.;
        deltaX[1] = plist[0]->at(i).x()-plist[0]->at(i2).x();
        deltaX[2] = plist[0]->at(i).y()-plist[0]->at(i2).y();
        deltaX[3] = plist[0]->at(i).z()-plist[0]->at(i2).z();
        double r = rInCMFrame(p1mu, p2mu, deltaX);

        double potential = CornellPotential(r);
          
        //Finally, the energy of the pair in the CM frame:
        double E_CM = M_inv+potential;

        //Hadronize the recombinant quarkonium:
        int status = 82;
        //cout << "Recombinant quarkonium should be formed, with status = " << status << endl;
        //cout << "tFinal of parton " << i << " = " << plist[0]->at(i).tFinal() << endl;
        //cout << "tFinal of parton " << antiI << " = " << plist[0]->at(antiI).tFinal() << endl;
        hadronizeByColorEvaporation(id1, id2, E_CM, p1mu, p2mu, M1, M2, status);  
        plist[0]->at(i2).frozen(1);
      }
        }
      }
    }

  return counter;
}

int MARTINI::hadronizeHeavyQuarks(vector<Parton> ** plist){
  pythia.event.reset(); // clear the event record to fill in the partons resulting from evolution below

  //If fullEvent == 0, this is relatively simple: we only determine the invariant mass of the pair 
  //in the CM frame and use this to hadronize the pairs into mesons:
  if(fullEvent == 0 || totalNNs == 1){
      int imax = plist[0]->size();
      for ( int i=0; i<imax; i++)                                   // loop over all partons
      {
        int id = plist[0]->at(i).id();
        if (abs(id) != 4 && abs(id) != 5) continue;

        if(plist[0]->at(i).frozen() != 0) continue;

        plist[0]->at(i).frozen(1);

        int antiI = plist[0]->at(i).antiI();

        if(plist[0]->at(antiI).frozen() == 0){
          int id1 = plist[0]->at(i).id();
          int id2 = plist[0]->at(antiI).id();

          //Any quarkonium produced will be the surviving component:
          int status = 81;

          //The relative momentum of the pair:
          Vec4 p1 = plist[0]->at(i).p();
          Vec4 p2 = plist[0]->at(antiI).p();

          double M1 = plist[0]->at(i).mass();
          double M2 = plist[0]->at(antiI).mass();

          double p1mu[4], p2mu[4];
          p1mu[1] = p1.px();
          p1mu[2] = p1.py();
          p1mu[3] = p1.pz();
          p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
          p2mu[1] = p2.px();
          p2mu[2] = p2.py();
          p2mu[3] = p2.pz();
          p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);

          double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
                  -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
                  -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
                  -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );

          double deltaX[4];
          deltaX[0] = 0.;
          deltaX[1] = plist[0]->at(i).x()-plist[0]->at(antiI).x();
          deltaX[2] = plist[0]->at(i).y()-plist[0]->at(antiI).y();
          deltaX[3] = plist[0]->at(i).z()-plist[0]->at(antiI).z();

          double r = rInCMFrame(p1mu, p2mu, deltaX);

          double potential = CornellPotential(r);

          //Finally, the energy of the pair in the CM frame:
          double E_CM = M_inv+potential;

          int bound = isBound(id1, id2, E_CM);

          if(bound == 1){
            hadronizeByColorEvaporation(id1, id2, E_CM, p1mu, p2mu, M1, M2, status);
            plist[0]->at(antiI).frozen(1);
            continue;
          }
        }

        Vec4 p1 = plist[0]->at(i).p();
        double M1 = plist[0]->at(i).mass();
        int id1 = plist[0]->at(i).id();
        double p1mu[4];
        p1mu[1] = p1.px();
        p1mu[2] = p1.py();
        p1mu[3] = p1.pz();
        p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
        hadronizeOneHeavyQuark(id1, p1mu, M1);
        continue;
      }

    //If the heavy quark did not evolve in the medium, fragment its momentum and hadronize it to 
    //an open heavy flavor meson:



//    int id1 = plist[0]->at(0).id();
//    int id2 = plist[0]->at(1).id();
//
//    //The status of surviving quarkonium is always 81:
//    int status = 81;
//
//    //The separation of the quarks in the lab frame:
//    double dx = plist[0]->at(0).x()-plist[0]->at(1).x();
//    double dy = plist[0]->at(0).y()-plist[0]->at(1).y();
//    double dz = plist[0]->at(0).z()-plist[0]->at(1).z();
//    double r = sqrt(dx*dx+dy*dy+dz*dz);
//    //After kinetic freezeout, the evolution of the heavy quarks is completely adiabatic. Therefore, 
//    //the Cornell potential should be used to determine the binding energy of the pair:
//    double potential;
//    potential = CornellPotential(r);
//
//    Vec4 p1 = plist[0]->at(0).p();
//    Vec4 p2 = plist[0]->at(1).p();
//    double M1 = plist[0]->at(0).mass();
//    double M2 = plist[0]->at(1).mass();
//
//    double p1mu[4], p2mu[4];
//    p1mu[1] = p1.px();
//    p1mu[2] = p1.py();
//    p1mu[3] = p1.pz();
//    p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
//    p2mu[1] = p2.px();
//    p2mu[2] = p2.py();
//    p2mu[3] = p2.pz();
//    p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);
//
//    //Let us just calculate with the invariant mass and the potential:
//    double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
//      -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
//      -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
//      -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );
//
//    double E_CM = M_inv+potential;
//
//    //The hadronization:
//    hadronizeByColorEvaporation(id1, id2, E_CM, p1mu, p2mu, M1, M2, status);
  }


//  if(fullEvent == 1 && totalNNs > 1 ){
//    //Determining recombinant production in a full event will be more difficult, but crucial. 
//    //First, the element of the permutation group describing the pairing of heavy quarks. 
//    //The quark from collision i, i=0, ..., totalNNs, is paired with the anti-quark in collision
//    //pairing[i]:
//    int pairing[totalNNs];
//    for(int ipair=0; ipair<totalNNs; ipair++){
//      pairing[ipair] = ipair;
//    }
//
//    //Whether or not the initial heavy quark pair is not hadronizing statistically. This can happen for two reasons:
//    //either the pair forms quarkonium and does not interact strongly with the medium, or the pair never 
//    //experienced in-medium evolution.
//    int notstatistical[totalNNs];
//    //int notevolved[totalNNs];
//    for(int ipair=0; ipair<totalNNs; ipair++){
//      notstatistical[totalNNs] = 0.;
//      //notevolved[totalNNs] = 0.;
//
//      //First, those pairs which do not hadronize statistically due to being bound:
//      int id1, id2;
//      id1 = plist[0]->at(2*ipair).id();
//      id2 = plist[0]->at(2*ipair+1).id();
//
//      double dx, dy, dz;
//      dx = plist[0]->at(2*ipair).x()-plist[0]->at(2*ipair+1).x();
//      dy = plist[0]->at(2*ipair).y()-plist[0]->at(2*ipair+1).y();
//      dz = plist[0]->at(2*ipair).z()-plist[0]->at(2*ipair+1).z();
//
//      double r = sqrt(dx*dx+dy*dy+dz*dz);
//      double potential = CornellPotential(r);
//
//      //The relative momentum of the pair:
//      Vec4 p1 = plist[0]->at(2*ipair).p();
//      Vec4 p2 = plist[0]->at(2*ipair+1).p();
//
//      double M1 = plist[0]->at(2*ipair).mass();
//      double M2 = plist[0]->at(2*ipair+1).mass();
//
//      double p1mu[4], p2mu[4];
//      p1mu[1] = p1.px();
//      p1mu[2] = p1.py();
//      p1mu[3] = p1.pz();
//      p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
//      p2mu[1] = p2.px();
//      p2mu[2] = p2.py();
//      p2mu[3] = p2.pz();
//      p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);
//
//      double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
//        -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
//        -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
//        -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );
//
//      //Finally, the energy of the pair in the CM frame:
//      double E_CM = M_inv+potential;
//
//      //All that for this:
//      if((abs(id1) ==4 && abs(id2) == 4) && E_CM<E_CCB_BOUND) notstatistical[ipair] = 1;
//      if((abs(id1) ==5 && abs(id2) == 5) && E_CM<E_BBB_BOUND) notstatistical[ipair] = 1;
//      if(( (abs(id1) ==4 && abs(id2) == 5) || (abs(id1) == 5 && abs(id2) == 4) ) && E_CM<E_BBC_BOUND) notstatistical[ipair] = 1;
//
//      //Finally, whether or not the pair evolved in medium:
//      //double z = 0.5*(plist[0]->at(2*ipair).z()+plist[0]->at(2*ipair+1).z());
//      //double t_final = 0.5*(plist[0]->at(2*ipair).tFinal()+plist[0]->at(2*ipair+1).tFinal() );
//      //double tau = sqrt(t_final*t_final-z*z);
//
//      //Only need to check the final time of one particle to see if both were evolved, when we are working with
//      //examineHQ==2:
//      double t_final = plist[0]->at(2*ipair).tFinal();
//      if(t_final == 0.){
//      notstatistical[ipair] = 1;
//      //notevolved[ipair] = 1;
//      }
//    }
//
//    //We perform the Metropolis algorithm to determine the pairing of these heavy quarks:
//    for(int iter=0; iter<N_SAMPLES; iter++){
//      //Select a random 2-cycle:
//      int cycle[2] ;
//      cycle[0] = gsl_rng_uniform_int(gsl_rand, totalNNs);
//      bool select_random = true;
//      while(select_random){
//        int random_int = gsl_rng_uniform_int(gsl_rand, totalNNs);
//        if(random_int != cycle[0]){
//          cycle[1] = random_int;
//          select_random = false;
//        }
//      }
//      
//      //Skip the rest of the loop if the pair is not hadronizing statistically:
//      if(notstatistical[cycle[0]]==1) continue;
//      if(notstatistical[cycle[1]]==1) continue;
//      
//      int pairing_new[totalNNs];
//      //The new pairing is mostly the same...
//      for(int ipair=0; ipair<totalNNs; ipair++){
//        pairing_new[ipair] = pairing[ipair];
//      }
//      //Except that it is multiplied on the left by the 2-cycle:
//      pairing_new[cycle[0]]=pairing[cycle[1]];
//      pairing_new[cycle[1]]=pairing[cycle[0]];
//      
//      //This is fine, but we need the exact locations of these particles in plist:
//      int i1, i2, antiI1, antiI2;
//      if(plist[0]->at(2*cycle[0]).id()>0) i1 = 2*cycle[0];
//      else i1 = 2*cycle[0]+1;
//      if(plist[0]->at(2*cycle[1]).id()>0) i2 = 2*cycle[1];
//      else i2 = 2*cycle[1]+1;
//      if(plist[0]->at(2*pairing[cycle[0]]).id()<0) antiI1 = 2*pairing[cycle[0]];
//      else antiI1 = 2*pairing[cycle[0]]+1;
//      if(plist[0]->at(2*pairing[cycle[1]]).id()<0) antiI2 = 2*pairing[cycle[1]];
//      else antiI2 = 2*pairing[cycle[1]]+1;
//      
//      //Determine whether or not to accept or reject this pairing,
//      //based on the Metropolis algorithm.
//      //First, determine the difference in the total potential energy for the new and old pairings:
//      double r1, r2, r1new, r2new;
//      double dV;
//      
//      //We are going to reuse these dx variables for calculating the various r's and ultimately, dV:
//      double dxtemp, dytemp, dztemp;
//      dxtemp = plist[0]->at(i1).x()-plist[0]->at(antiI1).x();
//      dytemp = plist[0]->at(i1).y()-plist[0]->at(antiI1).y();
//      dztemp = plist[0]->at(i1).z()-plist[0]->at(antiI1).z();
//      r1 = sqrt(dxtemp*dxtemp+dytemp*dytemp+dztemp*dztemp);
//      dxtemp = plist[0]->at(i2).x()-plist[0]->at(antiI2).x();
//      dytemp = plist[0]->at(i2).y()-plist[0]->at(antiI2).y();
//      dztemp = plist[0]->at(i2).z()-plist[0]->at(antiI2).z();
//      r2 = sqrt(dxtemp*dxtemp+dytemp*dytemp+dztemp*dztemp);
//      dxtemp = plist[0]->at(i1).x()-plist[0]->at(antiI2).x();
//      dytemp = plist[0]->at(i1).y()-plist[0]->at(antiI2).y();
//      dztemp = plist[0]->at(i1).z()-plist[0]->at(antiI2).z();
//      r1new = sqrt(dxtemp*dxtemp+dytemp*dytemp+dztemp*dztemp);
//      dxtemp = plist[0]->at(i2).x()-plist[0]->at(antiI1).x();
//      dytemp = plist[0]->at(i2).y()-plist[0]->at(antiI1).y();
//      dztemp = plist[0]->at(i2).z()-plist[0]->at(antiI1).z();
//      r2new = sqrt(dxtemp*dxtemp+dytemp*dytemp+dztemp*dztemp);
//
//      //cout << "r1 = " << r1 << ", r2 = " << r2 << ", r1new = " << r1new << ", r2new = " << r2new << endl;
//
//      dV = CornellPotential(r1new)+CornellPotential(r2new)-CornellPotential(r1)-CornellPotential(r2);
//      double BoltzmannFactor = exp(-dV/T_C_HQ);
//
//      //The Metropolis algorithm: if the Boltzmann factor is greater than unity, immediately accept the permutation
//      //(this is arbitrary but determines the next step):
//      //cout << "dV = " << dV << endl;
//      //cout << "T_C_HQ = " << T_C_HQ << endl;
//      //cout << "Boltzmann factor = " << BoltzmannFactor << endl;
//      if(BoltzmannFactor > 1){
//         for(int ipair=0; ipair<totalNNs; ipair++){
//           pairing[ipair] = pairing_new[ipair];
//         }
//      }
//      //If it is less than unity, accept the permutation at a rate that maintains detailed balance:
//      else{
//         double rdouble = gsl_rng_uniform(gsl_rand);
//         //cout << "rdouble = " << rdouble << endl;
//         if(rdouble < BoltzmannFactor){
//           for(int ipair=0; ipair<totalNNs; ipair++){
//             pairing[ipair] = pairing_new[ipair];
//           }
//         }
//      }
//      ////Output the new permutation for testing:
//      //for(int ipair=0; ipair<totalNNs; ipair++){
//      //cout << pairing[ipair] << " ";
//      //}
//      //cout << endl;
//    }
//
//    ////Output the new permutation for testing:
//    //cout << "The new permutation: " ;
//    //for(int ipair=0; ipair<totalNNs; ipair++){
//    //cout << pairing[ipair] << " ";
//    //}
//    //cout << endl;
//
//    //Now that the pairings have been sampled, we loop through the full event, adding J/Psi's, 
//    //Upsilons, b\bar{c} states, c\bar{b} states, and unhadronized free heavy quarks to the 
//    //event record:
//    for(int ie=0; ie<totalNNs; ie++){
//      int i1, i2;
//      if(plist[0]->at(2*ie).id()>0) i1 = 2*ie;
//      else i1 = 2*ie+1;
//      if(plist[0]->at(2*pairing[ie]).id()<0) i2 = 2*pairing[ie];
//      else i2 = 2*pairing[ie]+1;
//      int id1 = plist[0]->at(i1).id();
//      int id2 = plist[0]->at(i2).id();
//      //Distinguish between diagonal and recombinant production:
//      int status;
//      if(ie == pairing[ie]) status = 81;
//      else status = 82;
//      
//      //Determine again the binding energy in the CM frame:
//      //The separation of the quarks in the lab frame:
//      double dx = plist[0]->at(i1).x()-plist[0]->at(i2).x();
//      double dy = plist[0]->at(i1).y()-plist[0]->at(i2).y();
//      double dz = plist[0]->at(i1).z()-plist[0]->at(i2).z();
//      double r = sqrt(dx*dx+dy*dy+dz*dz);
//      //After kinetic freezeout, the evolution of the heavy quarks is completely adiabatic. Therefore, 
//      //the Cornell potential should be used to determine the binding energy of the pair:
//      double potential;
//      //if(notevolved[ie] == 1){
//      //potential = CornellPotential(0.);
//      //}
//      potential = CornellPotential(r);
//
//      //The relative momentum of the pair:
//      Vec4 p1 = plist[0]->at(i1).p();
//      Vec4 p2 = plist[0]->at(i2).p();
//      double M1 = plist[0]->at(i1).mass();
//      double M2 = plist[0]->at(i2).mass();
//
//      double p1mu[4], p2mu[4];
//      p1mu[1] = p1.px();
//      p1mu[2] = p1.py();
//      p1mu[3] = p1.pz();
//      p1mu[0] = sqrt(M1*M1+p1mu[1]*p1mu[1]+p1mu[2]*p1mu[2]+p1mu[3]*p1mu[3]);
//      p2mu[1] = p2.px();
//      p2mu[2] = p2.py();
//      p2mu[3] = p2.pz();
//      p2mu[0] = sqrt(M2*M2+p2mu[1]*p2mu[1]+p2mu[2]*p2mu[2]+p2mu[3]*p2mu[3]);
//
//      double M_inv = sqrt((p1mu[0]+p2mu[0])*(p1mu[0]+p2mu[0])
//        -(p1mu[1]+p2mu[1])*(p1mu[1]+p2mu[1])
//        -(p1mu[2]+p2mu[2])*(p1mu[2]+p2mu[2])
//        -(p1mu[3]+p2mu[3])*(p1mu[3]+p2mu[3]) );
//
//      //Finally, the energy of the pair in the CM frame:
//      double E_CM = M_inv+potential;
//
//      hadronizeByColorEvaporation(id1, id2, E_CM, p1mu, p2mu, M1, M2, status);
//    }
//  }
  return 0;
} 

// Uses the rejection method:
double MARTINI::SamplePetersonFunction(double M){
  double eps_Q;
  eps_Q = (CHARM_MASS/M)*(CHARM_MASS/M)*EPS_C;
  
  double z;
  bool zNotDetermined = true;
  while(zNotDetermined){

    double zSampled, fSampled;
    zSampled = gsl_rng_uniform(gsl_rand);
    fSampled = gsl_rng_uniform(gsl_rand)/eps_Q;

    double fValue = 1./(zSampled*(1.-1./zSampled-eps_Q/(1.-zSampled))*(1.-1./zSampled-eps_Q/(1.-zSampled)) );
    if(fSampled < fValue){
      z = zSampled;
      zNotDetermined = false;
    }
  }

  return z;
}

//Checks to see if a pair is bound:
int MARTINI::isBound(int id1, int id2, double E_CM){
  int isBound = 0;
  if(abs(id1)==4 && abs(id2)==4){
    if(E_CM < E_CCB_BOUND){
      isBound = 1;
    }
  }
  if(abs(id1)==5 && abs(id2)==5){
    if(E_CM < E_BBB_BOUND){
      isBound = 1;
    }
  }
  if((abs(id1)==4 && abs(id2)==5) || (abs(id1)==5 && abs(id2)==4) ){
    if(E_CM < E_BBC_BOUND){
      isBound = 1;
    }
  }
  return isBound;
}

int MARTINI::check_coalescence(int id1, double *pmu, double M1, double x, double y, double z, double t){
    double tau;

    HydroInfo hydroInfo;// antihydroInfo;
    PreHydroInfo prehydroInfo;// antihydroInfo;
  
    if(z*z<t*t) tau = sqrt(t*t-z*z); else return false;    // determine tau (z*z) should always be less than (t*t)
    double umu[4], pmuRest[4];
    if(tau<hydroTau0 && tau > prehydroTau0 && prehydroevolution == 1)
    {
      prehydroInfo = hydroSetup->getPreHydroValues(x, y, z, t);
      double eta_local = 0.5*log((t+z)/(t-z));
      double cosh_eta = cosh(eta_local);
      double sinh_eta = sinh(eta_local);
      umu[0] = prehydroInfo.utau*cosh_eta;  // gamma factor
      umu[1] = prehydroInfo.ux; 
      umu[2] = prehydroInfo.uy; 
      umu[3] = prehydroInfo.utau*sinh_eta;
      LorentzBooster(umu, pmu, pmuRest);                      // Boost into the fluid's rest frame
    } else if (tau >= hydroTau0)
    {
      hydroInfo = hydroSetup->getHydroValues(x, y, z, t); 
    
      double vx = hydroInfo.vx;
      double vy = hydroInfo.vy;
      double vz = hydroInfo.vz;
      double beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of flow velocity in units of c
      double gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
      umu[0] = gamma;
      umu[1] = gamma*vx;
      umu[2] = gamma*vy;
      umu[3] = gamma*vz;
      LorentzBooster(umu, pmu, pmuRest);                      // Boost into the fluid's rest frame
    }

    double pHQ = sqrt(pmuRest[1]*pmuRest[1] + pmuRest[2]*pmuRest[2] + pmuRest[3]*pmuRest[3]);
    double probability_to_coalesce = hqRates->get_coalesence_prob(pHQ);
    double probability_to_coalesce_meson = hqRates->get_coalesence_meson_prob(pHQ);
    double random = gsl_rng_uniform(gsl_rand);

    int decision = 0;

    if (random < coal_N_factor*probability_to_coalesce) decision = 2; //Coalescening
    if (random < coal_N_factor*probability_to_coalesce && random < coal_N_factor*probability_to_coalesce_meson) decision = 1;//Coalescening to meson

    return decision;
}

bool MARTINI::hadronize_by_coalescence(int id1, double *pHQ, double M, double x, double y, double z, double t){
   Vec4 pL;
   double pLight[4];
   bool coalesence_complete = false;
   int counter = 0;

   while(!coalesence_complete && counter < 100000)
   {
     counter++;
     pL = random->thermal_mass(T, 1, M);
     pLight[1] = pL.px();
     pLight[2] = pL.py();
     pLight[3] = pL.pz();
     pLight[0] = sqrt(M*M+pLight[1]*pLight[1]+pLight[2]*pLight[2]+pLight[3]*pLight[3]); //vecp loaded into an array for ease with LorentzBooster

     double Coal_normal = hqRates->get_Coal_normalization();

     double tau;

     HydroInfo hydroInfo;// antihydroInfo;
     PreHydroInfo prehydroInfo;// antihydroInfo;
  
     if(z*z<t*t) tau = sqrt(t*t-z*z); else break;    // determine tau (z*z) should always be less than (t*t)
     double umu[4], pLLab[4];
     if(tau<hydroTau0 && tau > prehydroTau0 && prehydroevolution == 1)
     {
       prehydroInfo = hydroSetup->getPreHydroValues(x, y, z, t);
       double eta_local = 0.5*log((t+z)/(t-z));
       double cosh_eta = cosh(eta_local);
       double sinh_eta = sinh(eta_local);
       umu[0] = prehydroInfo.utau*cosh_eta;  // gamma factor
       umu[1] = -1.*prehydroInfo.ux; 
       umu[2] = -1.*prehydroInfo.uy; 
       umu[3] = prehydroInfo.utau*sinh_eta;
       LorentzBooster(umu, pLight, pLLab);                      // Boost into the fluid's rest frame
     } else if (tau >= hydroTau0)
     {
       hydroInfo = hydroSetup->getHydroValues(x, y, z, t); 
     
       double vx = hydroInfo.vx;
       double vy = hydroInfo.vy;
       double vz = hydroInfo.vz;
       double beta = sqrt(vx*vx+vy*vy+vz*vz);                           // absolute value of flow velocity in units of c
       double gamma = 1./sqrt(1.-beta*beta);                            // gamma factor
       umu[0] = gamma;
       umu[1] = -1.*gamma*vx;
       umu[2] = -1.*gamma*vy;
       umu[3] = -1.*gamma*vz;
       LorentzBooster(umu, pLight, pLLab);                      // Boost into the fluid's rest frame
     }

     double v_cm[4];
     for (int ii = 1; ii <= 3; ii++)
     {
       v_cm[ii] = (pLLab[ii] + pHQ[ii])/(pLLab[0]+pHQ[0]);
     }
     double gamma = 1./(sqrt(1. - v_cm[1]*v_cm[1] - v_cm[2]*v_cm[2] - v_cm[3]*v_cm[3]));

     double pLQ_cm[4], pHQ_cm[4];

     for (int ii = 1; ii <= 3; ii++)
     {
       pLQ_cm[ii] = gamma*(pLLab[ii] - pLLab[0]*v_cm[ii]);
       pHQ_cm[ii] = gamma*(pHQ[ii]   - pHQ[0]*v_cm[ii]);
     }
     pLQ_cm[0] = gamma*(pLLab[0] - pLLab[1]*v_cm[1] - pLLab[2]*v_cm[2] - pLLab[3]*v_cm[3]);
     pHQ_cm[0] = gamma*(pHQ[0]   - pHQ[1]*v_cm[1]   - pHQ[2]*v_cm[2]   - pHQ[3]*v_cm[3]);

     double relK[4];
     for (int ii = 1; ii <= 3; ii++)
     {
       relK[ii] = (pHQ_cm[0]*pLQ_cm[ii] - pLQ_cm[0]*pHQ_cm[ii])/(pLQ_cm[0]+pHQ_cm[0]);
     }
     double relKsq = relK[1]*relK[1] + relK[2]*relK[2] + relK[3]*relK[3];

     double sigma = coal_sigma;

     double prob = exp(-relKsq*sigma*sigma);
     prob *= Coal_normal*coal_N_factor*pow(2.*sqrt(M_PI)*sigma, 3.)/36.;

     double random_num = gsl_rng_uniform(gsl_rand);

     if( random_num < prob)
     {
       coalesence_complete = true;
       double pnew[4];

       for (int ii = 0; ii < 4; ii++)
       {
         pnew[ii] = pHQ[ii] + pLLab[ii];
       }

       Vec4 pmeson;

       pmeson.p(pnew[1],pnew[2], pnew[3], pnew[0]);
       //For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
       //fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
       pythia.event.append(id1*411/abs(id1), 81, 0, 0, pmeson, DMASS);
     }

   }
   return coalesence_complete;
}

//void MARTINI::hadronizeOneHeavyQuark(int id, double *pmu, double M){
//
//  double z;
//  z = SamplePetersonFunction(M);
//
//  double pnew[4];
//  pnew[1] = z*pmu[1];
//  pnew[2] = z*pmu[2];
//  pnew[3] = z*pmu[3];
//
//  Vec4 pnewvec;
//
//  if(abs(id)==4){
//    pnew[0] = sqrt(pnew[1]*pnew[1]+pnew[2]*pnew[2]+pnew[3]*pnew[3]+DMASS*DMASS);
//    pnewvec.p(pnew[1],pnew[2], pnew[3], pnew[0]);
//    pythia.event.append(id*411/abs(id), 81, 0, 0, pnewvec, DMASS);
//  }
//  else{
//    pnew[0] = sqrt(pnew[1]*pnew[1]+pnew[2]*pnew[2]+pnew[3]*pnew[3]+BMASS*BMASS);
//    pnewvec.p(pnew[1],pnew[2], pnew[3], pnew[0]);
//    pythia.event.append(id*511/abs(id), 81, 0, 0, pnewvec, BMASS);
//  }
//}

void MARTINI::hadronizeOneHeavyQuark(int id, double *pmu, double M){

  double pmag = sqrt(pmu[1]*pmu[1] + pmu[2]*pmu[2] + pmu[3]*pmu[3]);

  // e1  = pT/|pT|, e2 = 
  double e1[3], e2[3], v[3];
  if (pmu[3]/pmag < 0.99) {
    v[0] = 0.; v[1] = 0.; v[2] = 1.;
  } else {
    v[0] = 0.; v[1] = 1.; v[2] = 0.;
  }
  e1[0] = (pmu[2]/pmag)*v[2] - (pmu[3]/pmag)*v[1];
  e1[1] = (pmu[3]/pmag)*v[0] - (pmu[1]/pmag)*v[2];
  e1[2] = (pmu[1]/pmag)*v[1] - (pmu[2]/pmag)*v[0];
  double e1_norm = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
  for (int ii = 0; ii < 3; ii++) e1[ii] /= e1_norm;

  e2[0] = (pmu[2]/pmag)*e1[2] - (pmu[3]/pmag)*e1[1];
  e2[1] = (pmu[3]/pmag)*e1[0] - (pmu[1]/pmag)*e1[2];
  e2[2] = (pmu[1]/pmag)*e1[1] - (pmu[2]/pmag)*e1[0];
  double e2_norm = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
  for (int ii = 0; ii < 3; ii++) e2[ii] /= e2_norm;

  double z;
  double sigma_kt = 0.335;//0.;//GeV
  z = SamplePetersonFunction(M);
  if(abs(id)==4){
    if (z*(pmu[0] + pmag) < DMASS) z = DMASS/(pmu[0] + pmag);
  } else if(abs(id)==5){
    if (z*(pmu[0] + pmag) < BMASS) z = BMASS/(pmu[0] + pmag);
  } else {
    return;
  }

  double p_plus_h = z*(pmu[0] + pmag);
  double kx = gsl_ran_gaussian(gsl_rand, sigma_kt);
  double ky = gsl_ran_gaussian(gsl_rand, sigma_kt);
  double kt_sq = kx*kx + ky*ky;
  double p_minus_h;
  //double p_minus_h = (DMASS*DMASS + kt_sq)/(2.*p_plus_h);
  if(abs(id)==4){
    p_minus_h = (DMASS*DMASS + kt_sq)/(p_plus_h);
  } else if(abs(id)==5){
    p_minus_h = (BMASS*BMASS + kt_sq)/(p_plus_h);
  } else {
    return;
  }

  double E_h_col  = 0.5*(p_plus_h + p_minus_h);
  double pz_h_col = 0.5*(p_plus_h - p_minus_h);

  double pnew[4];
  pnew[1] = (pmu[1]/pmag)*pz_h_col + kx*e1[0] + ky*e2[0];
  pnew[2] = (pmu[2]/pmag)*pz_h_col + kx*e1[1] + ky*e2[1];
  pnew[3] = (pmu[3]/pmag)*pz_h_col + kx*e1[2] + ky*e2[2];

  Vec4 pnewvec;

  if(abs(id)==4){
    pnew[0] = sqrt(pnew[1]*pnew[1]+pnew[2]*pnew[2]+pnew[3]*pnew[3]+DMASS*DMASS);
    //cout << pnew[1] << " " << z*pmu[1] << " " << pnew[2] << " " << z*pmu[2] << " " << pnew[3] << " " << z*pmu[3] << " " << pnew[0] << " " << sqrt(DMASS*DMASS+z*z*pmu[1]*pmu[1]+z*z*pmu[2]*pmu[2]+z*z*pmu[3]*pmu[3]) << endl;
    //if (sqrt(z*z*pmu[1]*pmu[1]+z*z*pmu[2]*pmu[2]) > 5. || sqrt(pnew[1]*pnew[1]+pnew[2]*pnew[2]) > 5.) cout << sqrt(pnew[1]*pnew[1]+pnew[2]*pnew[2]) << " " << sqrt(z*z*pmu[1]*pmu[1]+z*z*pmu[2]*pmu[2]) << endl;
    //cout << pmag << " " << e1[0] << " " << e1[1] << " " << e1[2] << "... " << e2[0] << " " << e2[1] << " " << e2[2] << "... " << pmu[1]/pmag << " " << pmu[2]/pmag << " " << pmu[3]/pmag << endl;
    pnewvec.p(pnew[1],pnew[2], pnew[3], pnew[0]);
    pythia.event.append(id*411/abs(id), 81, 0, 0, pnewvec, DMASS);
  }
  else if (abs(id) == 5) {
    pnew[0] = sqrt(pnew[1]*pnew[1]+pnew[2]*pnew[2]+pnew[3]*pnew[3]+BMASS*BMASS);
    pnewvec.p(pnew[1],pnew[2], pnew[3], pnew[0]);
    pythia.event.append(id*511/abs(id), 81, 0, 0, pnewvec, BMASS);
  }
}

//Given the necessary parameters, determines the quarkonium production with a (simplified) color evaporation model.
//This operates in two modes: in mode 1, the heavy quark with momentum p2mu hadronizes only if it is a constituent of quarkonium 
//with the heavy quark with momentum p1mu. In mode 2, both heavy quarks hadronize, whether quarkonium or open heavy flavor. For 
//this reason, the momenta are non-commutative!:
void MARTINI::hadronizeByColorEvaporation(int id1, int id2, double E_CM, double *p1mu, double *p2mu, double M1, double M2, int status){

  if(abs(id1)==4 && abs(id2)==4){
    //If the binding energy falls within certain ranges, either open charm or charmonium states 
    //are formed. We distinguish between three cases: pairs which we say hadronize into J/Psi particles
    //(and therefore have particle ID 443), particles which are bound yet which we will not specify at this 
    //point their state (which have the special particle ID ), and open heavy flavor mesons, whose momenta are 
    //determined with Peterson fragmentation:
    if(E_CM < E_JPSI){
      double pnewJPsi[4];
      pnewJPsi[1] = MJPsi*(p1mu[1]+p2mu[1])/E_CM;
      pnewJPsi[2] = MJPsi*(p1mu[2]+p2mu[2])/E_CM;
      pnewJPsi[3] = MJPsi*(p1mu[3]+p2mu[3])/E_CM;
      pnewJPsi[0] = sqrt(MJPsi*MJPsi+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
      Vec4 pnewvec;
      pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
      pythia.event.append(443, status, 0, 0, pnewvec, MJPsi);
    }
    if(E_CM > E_JPSI && E_CM < E_CCB_BOUND){
      double pnewJPsi[4];
      pnewJPsi[1] = MExcitedJPsi*(p1mu[1]+p2mu[1])/E_CM;
      pnewJPsi[2] = MExcitedJPsi*(p1mu[2]+p2mu[2])/E_CM;
      pnewJPsi[3] = MExcitedJPsi*(p1mu[3]+p2mu[3])/E_CM;
      pnewJPsi[0] = sqrt(MExcitedJPsi*MExcitedJPsi+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
      Vec4 pnewvec;
      pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
      pythia.event.append(443000, status, 0, 0, pnewvec, MExcitedJPsi);
    }
    if(E_CM > E_CCB_BOUND){
      double z1, z2;
      z1 = SamplePetersonFunction(M1);
      z2 = SamplePetersonFunction(M2);
      double pnew1[4], pnew2[4];
      pnew1[1] = z1*p1mu[1];
      pnew1[2] = z1*p1mu[2];
      pnew1[3] = z1*p1mu[3];
      pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+DMASS*DMASS);
      pnew2[1] = z2*p2mu[1];
      pnew2[2] = z2*p2mu[2];
      pnew2[3] = z2*p2mu[3];
      pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+DMASS*DMASS);
      Vec4 pnew1vec, pnew2vec;
      pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
      pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
      //For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
      //fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
      pythia.event.append(id1*411/abs(id1), 81, 0, 0, pnew1vec, DMASS);
      pythia.event.append(id2*411/abs(id2), 81, 0, 0, pnew2vec, DMASS);
    }
  }
  
  if(abs(id1)==5 && abs(id2)==5){
    if(E_CM < E_Y){
      double pnewJPsi[4];
      pnewJPsi[1] = MUpsilon*(p1mu[1]+p2mu[1])/E_CM;
      pnewJPsi[2] = MUpsilon*(p1mu[2]+p2mu[2])/E_CM;
      pnewJPsi[3] = MUpsilon*(p1mu[3]+p2mu[3])/E_CM;
      pnewJPsi[0] = sqrt(MUpsilon*MUpsilon+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
      Vec4 pnewvec;
      pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
      pythia.event.append(553, status, 0, 0, pnewvec, MUpsilon);
    }
    if(E_CM > E_Y && E_CM < E_BBB_BOUND){
      double pnewJPsi[4];
      pnewJPsi[1] = MExcitedUpsilon*(p1mu[1]+p2mu[1])/E_CM;
      pnewJPsi[2] = MExcitedUpsilon*(p1mu[2]+p2mu[2])/E_CM;
      pnewJPsi[3] = MExcitedUpsilon*(p1mu[3]+p2mu[3])/E_CM;
      pnewJPsi[0] = sqrt(MExcitedUpsilon*MExcitedUpsilon+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
      Vec4 pnewvec;
      pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
      pythia.event.append(553000, status, 0, 0, pnewvec, MExcitedUpsilon);
    }
    if(E_CM > E_BBB_BOUND){
      double z1, z2;
      z1 = SamplePetersonFunction(M1);
      z2 = SamplePetersonFunction(M2);
      //cout << "z1 = " << z1 << ", z2 = " << z2 << endl;
      double pnew1[4], pnew2[4];
      pnew1[1] = z1*p1mu[1];
      pnew1[2] = z1*p1mu[2];
      pnew1[3] = z1*p1mu[3];
      pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+BMASS*BMASS);
      pnew2[1] = z2*p2mu[1];
      pnew2[2] = z2*p2mu[2];
      pnew2[3] = z2*p2mu[3];
      pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+BMASS*BMASS);
      Vec4 pnew1vec, pnew2vec;
      pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
      pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
      //For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
      //fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
      pythia.event.append(id1*511/abs(id1), 81, 0, 0, pnew1vec, BMASS);
      pythia.event.append(id2*511/abs(id2), 81, 0, 0, pnew2vec, BMASS);
    }
  }
  if((abs(id1)==4 && abs(id2)==5) || (abs(id1)==5 && abs(id2)==4) ){
    if(E_CM < E_BBC){
      double pnewJPsi[4];
      pnewJPsi[1] = MB_c*(p1mu[1]+p2mu[1])/E_CM;
      pnewJPsi[2] = MB_c*(p1mu[2]+p2mu[2])/E_CM;
      pnewJPsi[3] = MB_c*(p1mu[3]+p2mu[3])/E_CM;
      pnewJPsi[0] = sqrt(MB_c*MB_c+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
      Vec4 pnewvec;
      pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
      pythia.event.append(543, status, 0, 0, pnewvec, MB_c);
    }
    if(E_CM > E_BBC && E_CM < E_BBC_BOUND){
      double pnewJPsi[4];
      pnewJPsi[1] = MExcitedB_c*(p1mu[1]+p2mu[1])/E_CM;
      pnewJPsi[2] = MExcitedB_c*(p1mu[2]+p2mu[2])/E_CM;
      pnewJPsi[3] = MExcitedB_c*(p1mu[3]+p2mu[3])/E_CM;
      pnewJPsi[0] = sqrt(MExcitedB_c*MExcitedB_c+pnewJPsi[1]*pnewJPsi[1]+pnewJPsi[2]*pnewJPsi[2]+pnewJPsi[3]*pnewJPsi[3]);
      Vec4 pnewvec;
      pnewvec.p(pnewJPsi[1],pnewJPsi[2],pnewJPsi[3],pnewJPsi[0]);
      pythia.event.append(543000, status, 0, 0, pnewvec, MExcitedB_c);
    }
    if(E_CM > E_BBC_BOUND){
      double z1, z2;
      z1 = SamplePetersonFunction(M1);
      z2 = SamplePetersonFunction(M2);  
      //cout << "z1 = " << z1 << ", z2 = " << z2 << endl;
      double pnew1[4], pnew2[4];
      
      pnew1[1] = z1*p1mu[1];
      pnew1[2] = z1*p1mu[2];
      pnew1[3] = z1*p1mu[3];
      
      pnew2[1] = z2*p2mu[1];
      pnew2[2] = z2*p2mu[2];
      pnew2[3] = z2*p2mu[3];
      
      //For now, we put all D-mesons into the pseudoscalar ground state; it may make sense at some point to 
      //fragment into D^+, D^-, D^*... based on the energy considerations of Peterson fragmentation:
      if(abs(id1)==4){
  pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+DMASS*DMASS);
  Vec4 pnew1vec;
  pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
  pythia.event.append(id1*411/abs(id1), 81, 0, 0, pnew1vec, DMASS);
      }
      else{
  pnew1[0] = sqrt(pnew1[1]*pnew1[1]+pnew1[2]*pnew1[2]+pnew1[3]*pnew1[3]+BMASS*BMASS);
  Vec4 pnew1vec;
  pnew1vec.p(pnew1[1],pnew1[2], pnew1[3], pnew1[0]);
  pythia.event.append(id1*511/abs(id1), 81, 0, 0, pnew1vec, BMASS);
      }
      if(abs(id2)==4){
  pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+DMASS*DMASS);
  Vec4 pnew2vec;
  pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
  pythia.event.append(id2*411/abs(id2), 81, 0, 0, pnew2vec, DMASS);
      }
      else{
  pnew2[0] = sqrt(pnew2[1]*pnew2[1]+pnew2[2]*pnew2[2]+pnew2[3]*pnew2[3]+BMASS*BMASS);
  Vec4 pnew2vec;
  pnew2vec.p(pnew2[1],pnew2[2], pnew2[3], pnew2[0]);
  pythia.event.append(id2*511/abs(id2), 81, 0, 0, pnew2vec, BMASS);
      }
    }
  }
  return;
}

//End rewrote routines for HQ by CFY -- MS
