#include "MARTINI.h"
#include <sstream>
#define numJets 100
#define nybins 200
#define DNDY2_MAX 10.
#define DNDY_MAX 10.
#define nptbins 50
#define nphibins 60
#define Y_MAX 1.0
#define Y_MAX_ALICE_v2 0.8
#define Y_MAX_dNdphi 2.0
#define Y_MAX_ALICE_RAA 0.5
#define Y_MAX_CMS 1.
#define PT_MAX 15.
#define R_MAX 10.
#define outputting true
#define nrapbins 200
#define ntau 200
//#define largeptbins 8

using namespace std;

int main(int argc, char* argv[]) 
{
  //The first bit of testing for heavy quarks!
  MARTINI martini;
  Event fullEvent;

  string setupFile = (argc > 1) ? argv[1] : "setup_main_JPsi_v2.dat";

  martini.readFile(setupFile);
  if (argc > 1)
  {
    int job_num = stoi(argv[1]);
    martini.init(job_num);
  } else {
    martini.init(0);
  }
  martini.settings.listAll();

  vector<Parton> ** plist;          // pointer to array of vector<Parton> objects
  vector<Parton> * plistInitial; 
  plist = new vector<Parton> *[2];  // pointer to array of vector<Parton> objects
  plistInitial = new vector<Parton>;
  plist[0] = new vector<Parton>;
  plist[1] = new vector<Parton>;

  string directory_name = "results/";
  string HQ_name = "LHC";

  int numEvents;
  int event_counter = 0;            // Counter of events

  int mt;                           // maximal time steps
  double maxTime=martini.returnMaxTime();
  double dtfm=martini.returnDtfm();
  mt = static_cast<int>(maxTime/dtfm+0.0001); // compute number of steps
  double zmax = maxTime;

  Vec4 itsP;                          // Its P

  //Define the histograms:
  double NJPsi_t[mt/10];
  double NJPsi_t_rec[mt/10];
  cout << "mt/10 = " << mt/10 << endl;

  double hydro_T[1000];
  double prehydro_T[100];
  int hydro_counter[1000];
  int prehydro_counter[100];
  
  for (int ind = 0; ind < 100; ind++){
    prehydro_T[ind] = 0.;
    prehydro_counter[ind] = 0;
  }

  for (int ind = 0; ind < 1000; ind++){
    hydro_T[ind] = 0.;
    hydro_counter[ind] = 0;
  }

//  double delpt[largeptbins];
//  double midpt[largeptbins];
//  if (PT_MAX > 100.) {
//	  cout << "Error:: Cannot track largeptbins for diffusion stability such a large ptmax" << endl;
//	  exit(1);
//  }
//  delpt[0] = 100. - PT_MAX;
//  delpt[1] = 900.;
//  delpt[2] = 9000.;
//  delpt[3] = 90000.;
//  delpt[4] = 900000.;
//  delpt[5] = 9000000.;
//  delpt[6] = 90000000.;
//  delpt[7] = 900000000.;
//  midpt[0] = (100. + PT_MAX)/2.;
//  midpt[1] = 550.;
//  midpt[2] = 5500.;
//  midpt[3] = 55000.;
//  midpt[4] = 550000.;
//  midpt[5] = 5500000.;
//  midpt[6] = 55000000.;
//  midpt[7] = 550000000.;
//  if (largeptbins > 8) {
//	  cout << "Error:: Cannot track largeptbins for diffusion stability such a large largeptbins" << endl;
//	  exit(1);
//  }
//
  double dN2piptdptdy_ccbar_beforeEvolution_CMS[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution_CMS[nptbins];
  double dN2piptdptdy_ccbar_beforeEvolution_ALICE_v2[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution_ALICE_v2[nptbins];
  double dN2piptdptdy_ccbar_beforeEvolution_ALICE_RAA[nptbins];
  double dN2piptdptdy_bbbar_beforeEvolution_ALICE_RAA[nptbins];
//  double dN2piptdptdy_ccbar_beforeEvolution_largept[largeptbins];
  
  double dN2piptdptdy_ccbar_CMS[nptbins];
  double dN2piptdptdy_bbbar_CMS[nptbins];
  double dN2piptdptdy_ccbar_ALICE_v2[nptbins];
  double dN2piptdptdy_bbbar_ALICE_v2[nptbins];
  double dN2piptdptdy_ccbar_ALICE_RAA[nptbins];
  double dN2piptdptdy_bbbar_ALICE_RAA[nptbins];
//  double dN2piptdptdy_ccbar_largept[largeptbins];
  double cos2Phi_dN2piptdptdy_ccbar_CMS[nptbins];
  double cos2Phi_dN2piptdptdy_bbbar_CMS[nptbins];
  double sin2Phi_dN2piptdptdy_ccbar_CMS[nptbins];
  double sin2Phi_dN2piptdptdy_bbbar_CMS[nptbins];
  double cos2Phi_dN2piptdptdy_ccbar_ALICE_v2[nptbins];
  double cos2Phi_dN2piptdptdy_bbbar_ALICE_v2[nptbins];
  double sin2Phi_dN2piptdptdy_ccbar_ALICE_v2[nptbins];
  double sin2Phi_dN2piptdptdy_bbbar_ALICE_v2[nptbins];
  double cos2Phi_dN2piptdptdy_ccbar_ALICE_RAA[nptbins];
  double cos2Phi_dN2piptdptdy_bbbar_ALICE_RAA[nptbins];
  double sin2Phi_dN2piptdptdy_ccbar_ALICE_RAA[nptbins];
  double sin2Phi_dN2piptdptdy_bbbar_ALICE_RAA[nptbins];

  double dN2piptdptdy_DDbar_CMS[nptbins];
  double dN2piptdptdy_DDbar_ALICE_v2[nptbins];
  double dN2piptdptdy_DDbar_ALICE_RAA[nptbins];
  double dN2piptdptdy_DBDBbar[nptbins];
  double cos2Phi_dN2piptdptdy_DDbar_CMS[nptbins];
  double cos2Phi_dN2piptdptdy_DDbar_ALICE_v2[nptbins];
  double cos2Phi_dN2piptdptdy_DDbar_ALICE_RAA[nptbins];
  double cos2Phi_dN2piptdptdy_DBDBbar[nptbins];
  double sin2Phi_dN2piptdptdy_DDbar_CMS[nptbins];
  double sin2Phi_dN2piptdptdy_DDbar_ALICE_v2[nptbins];
  double sin2Phi_dN2piptdptdy_DDbar_ALICE_RAA[nptbins];
  double sin2Phi_dN2piptdptdy_DBDBbar[nptbins];

  double dN2piptdptdy_JPsi[nptbins];
  double dN2piptdptdy_ExcitedJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_JPsi[nptbins];
  double cos2Phi_dN2piptdptdy_ExcitedJPsi[nptbins];
  double sin2Phi_dN2piptdptdy_JPsi[nptbins];
  double sin2Phi_dN2piptdptdy_ExcitedJPsi[nptbins];
  double dN2piptdptdy_OnlyDirectJPsi[nptbins];
  double dN2piptdptdy_OnlyDirectExcitedJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_OnlyDirectJPsi[nptbins];
  double cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[nptbins];
  double sin2Phi_dN2piptdptdy_OnlyDirectJPsi[nptbins];
  double sin2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[nptbins];

  double dN2piptdptdy_Upsilon[nptbins];
  double dN2piptdptdy_ExcitedUpsilon[nptbins];

  double dN2piptdptdy_bcbar[nptbins];
  double dN2piptdptdy_Excitedbcbar[nptbins];

  //Counting the number of bound states:
  double NJPsi, NExcitedJPsi, NUpsilon, NExcitedUpsilon, NBc, NExcitedBc;
  NJPsi = NExcitedJPsi = 0.;
  NUpsilon = NExcitedUpsilon = 0.;
  NBc = NExcitedBc = 0.;

//  for(int ip=0; ip<largeptbins; ip++){
//    dN2piptdptdy_ccbar_beforeEvolution_largept[ip] = 0.;
//    dN2piptdptdy_ccbar_beforeHydro_largept[ip] = 0.;
//    dN2piptdptdy_ccbar_largept[ip] = 0.;
//  }

  for(int ip=0; ip<nptbins; ip++){
    dN2piptdptdy_ccbar_beforeEvolution_CMS[ip] = 0.;
    dN2piptdptdy_bbbar_beforeEvolution_CMS[ip] = 0.;

    dN2piptdptdy_ccbar_beforeEvolution_ALICE_v2[ip] = 0.;
    dN2piptdptdy_bbbar_beforeEvolution_ALICE_v2[ip] = 0.;

    dN2piptdptdy_ccbar_beforeEvolution_ALICE_RAA[ip] = 0.;
    dN2piptdptdy_bbbar_beforeEvolution_ALICE_RAA[ip] = 0.;

//    dN2piptdptdy_ccbar_beforeHydro[ip] = 0.;
//    dN2piptdptdy_bbbar_beforeHydro[ip] = 0.;

    dN2piptdptdy_ccbar_CMS[ip] = 0.;
    dN2piptdptdy_bbbar_CMS[ip] = 0.;
    cos2Phi_dN2piptdptdy_ccbar_CMS[ip] = 0.;
    cos2Phi_dN2piptdptdy_bbbar_CMS[ip] = 0.;
    sin2Phi_dN2piptdptdy_ccbar_CMS[ip] = 0.;
    sin2Phi_dN2piptdptdy_bbbar_CMS[ip] = 0.;

    dN2piptdptdy_ccbar_ALICE_v2[ip] = 0.;
    dN2piptdptdy_bbbar_ALICE_v2[ip] = 0.;
    cos2Phi_dN2piptdptdy_ccbar_ALICE_v2[ip] = 0.;
    cos2Phi_dN2piptdptdy_bbbar_ALICE_v2[ip] = 0.;
    sin2Phi_dN2piptdptdy_ccbar_ALICE_v2[ip] = 0.;
    sin2Phi_dN2piptdptdy_bbbar_ALICE_v2[ip] = 0.;

    dN2piptdptdy_ccbar_ALICE_RAA[ip] = 0.;
    dN2piptdptdy_bbbar_ALICE_RAA[ip] = 0.;
    cos2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ip] = 0.;
    cos2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ip] = 0.;
    sin2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ip] = 0.;
    sin2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ip] = 0.;

    dN2piptdptdy_DDbar_CMS[ip] = 0.;
    cos2Phi_dN2piptdptdy_DDbar_CMS[ip] = 0.;
    sin2Phi_dN2piptdptdy_DDbar_CMS[ip] = 0.;

    dN2piptdptdy_DDbar_ALICE_v2[ip] = 0.;
    cos2Phi_dN2piptdptdy_DDbar_ALICE_v2[ip] = 0.;
    sin2Phi_dN2piptdptdy_DDbar_ALICE_v2[ip] = 0.;

    dN2piptdptdy_DDbar_ALICE_RAA[ip] = 0.;
    cos2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ip] = 0.;
    sin2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ip] = 0.;

    dN2piptdptdy_DBDBbar[ip] = 0.;
    cos2Phi_dN2piptdptdy_DBDBbar[ip] = 0.;
    sin2Phi_dN2piptdptdy_DBDBbar[ip] = 0.;

    dN2piptdptdy_JPsi[ip] = 0.;
    dN2piptdptdy_ExcitedJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_JPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_ExcitedJPsi[ip] = 0.;
    sin2Phi_dN2piptdptdy_JPsi[ip] = 0.;
    sin2Phi_dN2piptdptdy_ExcitedJPsi[ip] = 0.;
    dN2piptdptdy_OnlyDirectJPsi[ip] = 0.;
    dN2piptdptdy_OnlyDirectExcitedJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ip] = 0.;
    cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip] = 0.;
    sin2Phi_dN2piptdptdy_OnlyDirectJPsi[ip] = 0.;
    sin2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip] = 0.;

    dN2piptdptdy_Upsilon[ip] = 0.;
    dN2piptdptdy_ExcitedUpsilon[ip] = 0.;

    dN2piptdptdy_bcbar[ip] = 0.;
    dN2piptdptdy_Excitedbcbar[ip] = 0.;
  }

  double dNdy_ccbar_beforeEvolution[nybins];
  double dNdy_ccbar[nybins];
  double dNdpseudorap_ccbar_beforeEvolution[nybins];
  double dNdpseudorap_ccbar[nybins];
  double dNdstrap_ccbar[nybins];
  double dNdy_DDbar[nybins];
  double dNdy_JPsi[nybins];
  double dNdy_ExcitedJPsi[nybins];
  double delpt2_preeq[nybins];
  double delpt2[nybins];
  double delpt2_nopz[nybins];
  double delpt2_nopz_preeq[nybins];
  double dely2[nybins];

  for(int iy=0; iy<nybins; iy++){
    dNdy_ccbar_beforeEvolution[iy] = 0.;
    dNdy_ccbar[iy] = 0.;
    dNdpseudorap_ccbar_beforeEvolution[iy] = 0.;
    dNdpseudorap_ccbar[iy] = 0.;
    dNdstrap_ccbar[iy] = 0.;
    dNdy_DDbar[iy] = 0.;
    dNdy_JPsi[iy] = 0.;
    dNdy_ExcitedJPsi[iy] = 0.;
    delpt2_preeq[iy] = 0.;
    delpt2[iy] = 0.;
    delpt2_nopz_preeq[iy] = 0.;
    delpt2_nopz[iy] = 0.;
    dely2[iy] = 0.;
  }

  double dN2dphidy_ccbar_beforeEvolution[nphibins];
  double dN2dphidy_ccbar[nphibins];
  double dN2dphidy_DDbar[nphibins];
//  double dN2dphidy_ccbar_mixed_beforeEvolution_ALICE_v2[nphibins];
//  double dN2dphidy_ccbar_mixed_ALICE_v2[nphibins];
//  double dN2dphidy_DDbar_mixed_ALICE_v2[nphibins];
  double ccbar_pairs_beforeevolution = 0.;
  double ccbar_pairs = 0.;
  double DDbar_pairs = 0.;
//  double ccbar_mixed_pairs_beforeevolution = 0.;
//  double ccbar_mixed_pairs = 0.;
//  double DDbar_mixed_pairs = 0.;

  for(int iphi=0; iphi<nphibins; iphi++){
    dN2dphidy_ccbar_beforeEvolution[iphi] = 0.;
    dN2dphidy_ccbar[iphi] = 0.;
    dN2dphidy_DDbar[iphi] = 0.;
    //dN2dphidy_ccbar_mixed_beforeEvolution_ALICE_v2[iphi] = 0.;
    //dN2dphidy_ccbar_mixed_ALICE_v2[iphi] = 0.;
    //dN2dphidy_DDbar_mixed_ALICE_v2[iphi] = 0.;
  }

  //double dN4pir2dr_ccbar_beforeEvolution[nptbins];
  //double dN4pir2dr_ccbar[nptbins];
  //for(int ir=0; ir<nptbins; ir++){
  //dN4pir2dr_ccbar_beforeEvolution[ir] = 0.;
  //dN4pir2dr_ccbar[ir] = 0.;
  //}

  ////Now, we test the fragmentation of the heavy quarks with 
  ////generateTestEvent_HQ:
  //ofstream NboundOut;
  //NboundOut.open("Nbound.dat");
  //for(int il=0; il<50; il++){
  //double Nbound = 0.;
  //int NSAMP = 10000;
  //double L = 0.7+0.5*(double)il/50.;
  //cout << "L = " << L << endl;
  //for(int isamp=0; isamp<NSAMP; isamp++){
  //  plist[0]->clear();
  //  plist[1]->clear();
  //  martini.generateTestEvent_HQ(plist[0], L);
  //  martini.fragmentation_HQ(plist);
  //  int eventSize = martini.pythia.event.size();
  //  //cout << "Event size in the fragmentation testing = " << eventSize << endl;
    //  for(int ie=0; ie<eventSize; ie++){
  //
  ////cout << "martini.pythia.event[" << ie << "].id() = " << martini.pythia.event[ie].id() << endl;
  //if(martini.pythia.event[ie].id() == 443){
  //  Nbound += 1./(double)NSAMP;
  //}
  //if(martini.pythia.event[ie].id() == 443000){
  //  Nbound += 1./(double)NSAMP;
  //}
  //  }
  //}
  //NboundOut << L << " " << Nbound << endl;
  //}
  //NboundOut.close();

  cout << numJets << endl;

  while(event_counter<numJets)      // The main loop
    {
      event_counter++;

      plist[0]->clear();          // clear the parton list
      plist[1]->clear();          // clear the parton list

      cout << "event_counter  " << event_counter << endl;
      martini.generateEventHeavyQuarks_w_colllist_OR_IPGfile(plist[0]); 
//      martini.generateEventHeavyQuarks(plist[0]); 

      int eventSize = plist[0]->size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	//rapidity is more appropriate here:
	double pseudorap = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double y         = 0.5*log( (sqrt(pt*pt+pz*pz+itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz+itsM*itsM) - pz) );

	//if (abs(plist[0]->at(ie).id()) == 4) cout << plist[0]->at(ie).id() << " " << itsP.px() << " " << itsP.py() << " " << itsP.pz() << " " << itsM << " " << y << endl;
	if (plist[0]->at(ie).id()==4)
        {
	  for (int je=0; je<eventSize; je++){
            if (plist[0]->at(je).id() != -4) continue;
	    Vec4 itsPj = plist[0]->at(je).p();
	    double ptj = sqrt(itsPj.px()*itsPj.px() + itsPj.py()*itsPj.py());
	    double itsMj = plist[0]->at(je).mass();
            double pzj = itsPj.pz();
            //rapidity is more appropriate here:
            double yj         = 0.5*log( (sqrt(ptj*ptj+pzj*pzj+itsMj*itsMj) + pzj)/(sqrt(ptj*ptj+pzj*pzj+itsMj*itsMj) - pzj) );
	    if(fabs(y)<Y_MAX_dNdphi && fabs(yj)<Y_MAX_dNdphi){
	      double phi1 = atan(itsP.py()/itsP.px());
	      if (itsP.px() < 0) phi1 = phi1 + M_PI;
	      double phi2 = atan(itsPj.py()/itsPj.px());
	      if (itsPj.px() < 0) phi2 = phi2 + M_PI;
	      double delphi = phi1-phi2;
	      if (delphi < 0) delphi = - delphi;
	      if (delphi > M_PI ) delphi = 2.*M_PI - delphi;
	      int iphi =(int)( delphi*(double)nphibins/(M_PI));
	      if (iphi == nphibins) iphi = nphibins - 1;
	      //if (iphi < 0 || iphi >= nphibins) std::cout << iphi << endl;
              dN2dphidy_ccbar_beforeEvolution[iphi] += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
	      ccbar_pairs_beforeevolution += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
	      //cout << "Same  " << phi1 << " " << phi2 << " " << delphi << endl;
            }
	    //if(fabs(y)<Y_MAX_ALICE_v2 && yj<(Y_MAX_ALICE_v2+3.) && yj>(3.-Y_MAX_ALICE_v2)){
	    //  double phi1 = atan(itsP.py()/itsP.px());
	    //  if (itsP.px() < 0) phi1 = phi1 + M_PI;
	    //  double phi2 = atan(itsPj.py()/itsPj.px());
	    //  if (itsPj.px() < 0) phi2 = phi2 + M_PI;
	    //  double delphi = phi1-phi2;
	    //  if (delphi < 0) delphi = - delphi;
	    //  if (delphi > M_PI ) delphi = 2.*M_PI - delphi;
	    //  int iphi =(int)( delphi*(double)nphibins/(M_PI));
	    //  if (iphi == nphibins) iphi = nphibins - 1;
	    //  //if (iphi < 0 || iphi >= nphibins) std::cout << iphi << endl;
            //  dN2dphidy_ccbar_mixed_beforeEvolution_ALICE_v2[iphi] += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
	    //  ccbar_mixed_pairs_beforeevolution += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
	    //  //cout << "Different  " << phi1 << " " << phi2 << " " << delphi << endl;
            //}
	  }
	}

	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy>=0 && iy<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdy_ccbar_beforeEvolution[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	}
	int ipseudorap = (int)( (pseudorap+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(ipseudorap>=0 && ipseudorap<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdpseudorap_ccbar_beforeEvolution[ipseudorap] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX_CMS && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_beforeEvolution_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_beforeEvolution_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }
	if(fabs(y)<Y_MAX_ALICE_v2 && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_beforeEvolution_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_beforeEvolution_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }
	if(fabs(y)<Y_MAX_ALICE_RAA && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_beforeEvolution_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_beforeEvolution_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }

      }

      double hydrotau0 = martini.returnHydroTau0();
      bool wrote_prehydro_list = false;

      int counter = 0;
      if (martini.returnEvolution() == 1)
	{
	  for(int i=0; i<mt; i++)
	    {
	      //counter = martini.evolveHeavyQuarks(plist, counter, i);
	      counter = martini.evolveAndHadronizeHeavyQuarks(plist, hydro_T, prehydro_T, hydro_counter, prehydro_counter, counter, i);
	      counter+=1;
	//      if (!wrote_prehydro_list && time >= hydrotau0)
	//	{
	//	  wrote_prehydro_list = true;
        //          int eventSize = plist[0]->size();

        //          for(int ie=0; ie<eventSize; ie++){
        //            Vec4 itsP = plist[0]->at(ie).p();
        //            double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
        //            double itsM = plist[0]->at(ie).mass();
        //            double pz = itsP.pz();
        //            //Pseudorapidity is more appropriate here:
        //            double y = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );

        //            //The dN/dy histograms:
        //            int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
        //            if(iy>=0 && iy<nybins){
        //              if(abs(martini.pythia.event[ie].id())==4){
        //                dNdy_ccbar_beforeHydro[iy] += (double)nybins/(2.*DNDY_MAX);
        //              }
        //            }

        //            if(fabs(y)<Y_MAX && pt<PT_MAX){
        //              int ipt = (int)((pt/PT_MAX)*nptbins );

        //              if(abs(martini.pythia.event[ie].id())==4){
        //                dN2piptdptdy_ccbar_beforeHydro[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
        //              }
        //              if(abs(martini.pythia.event[ie].id())==5){
        //                dN2piptdptdy_bbbar_beforeHydro[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
        //              }
	//            } else if (fabs(y)<Y_MAX) {
	//              int ipt = (int)(log(pt)/log(10.)) - 1;
	//              if(abs(martini.pythia.event[ie].id())==4 && ipt < largeptbins){
	//                dN2piptdptdy_ccbar_beforeHydro_largept[ipt] += delpt[ipt]*(0.5/Y_MAX)/(2.*M_PI*midpt[ipt] );
	//              }
        //            }

        //          }

	//	}
	    }
	}

      eventSize = plist[0]->size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = plist[0]->at(ie).p();
	Vec4 itsIP = plist[0]->at(ie).Initialp();
	Vec4 itsEOPP = plist[0]->at(ie).EndOfPreeqp();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = plist[0]->at(ie).mass();
	double pz = itsP.pz();
	//rapidity is more appropriate here:
	double pseudorap = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double y = 0.5*log( (sqrt(pt*pt+pz*pz+itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz+itsM*itsM) - pz) );
	double cos2Phi = 0.;
	double sin2Phi = 0.;
	if(pt > 1e-5){
	  cos2Phi = (itsP.px()*itsP.px() - itsP.py()*itsP.py())/(pt*pt);
	  sin2Phi = 2.*(itsP.px()*itsP.py())/(pt*pt);
	}

	if (plist[0]->at(ie).id()==4)
        {
	  for (int je=0; je<eventSize; je++){
            if (plist[0]->at(je).id() != -4) continue;
	    Vec4 itsPj = plist[0]->at(je).p();
	    double ptj = sqrt(itsPj.px()*itsPj.px() + itsPj.py()*itsPj.py());
	    double itsMj = plist[0]->at(je).mass();
            double pzj = itsPj.pz();
            //rapidity is more appropriate here:
            double yj         = 0.5*log( (sqrt(ptj*ptj+pzj*pzj+itsMj*itsMj) + pzj)/(sqrt(ptj*ptj+pzj*pzj+itsMj*itsMj) - pzj) );
	    if(fabs(y)<Y_MAX_dNdphi && fabs(yj)<Y_MAX_dNdphi){
	      double phi1 = atan(itsP.py()/itsP.px());
	      if (itsP.px() < 0) phi1 = phi1 + M_PI;
	      double phi2 = atan(itsPj.py()/itsPj.px());
	      if (itsPj.px() < 0) phi2 = phi2 + M_PI;
	      double delphi = phi1-phi2;
	      if (delphi < 0) delphi = - delphi;
	      if (delphi > M_PI) delphi = 2.*M_PI - delphi;
	      int iphi =(int)( delphi*(double)nphibins/(M_PI));
	      if (iphi == nphibins) iphi = nphibins - 1;
	      //if (iphi < 0 || iphi >= nphibins) std::cout << iphi << endl;
              dN2dphidy_ccbar[iphi] += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
	      ccbar_pairs += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
            }
	    //if(fabs(y)<Y_MAX_ALICE_v2 && yj<(Y_MAX_ALICE_v2+3.) && yj>(3.-Y_MAX_ALICE_v2)){
	    //  double phi1 = atan(itsP.py()/itsP.px());
	    //  if (itsP.px() < 0) phi1 = phi1 + M_PI;
	    //  double phi2 = atan(itsPj.py()/itsPj.px());
	    //  if (itsPj.px() < 0) phi2 = phi2 + M_PI;
	    //  double delphi = phi1-phi2;
	    //  if (delphi < 0) delphi = - delphi;
	    //  if (delphi > M_PI) delphi = 2.*M_PI - delphi;
	    //  int iphi =(int)( delphi*(double)nphibins/(M_PI));
	    //  if (iphi == nphibins) iphi = nphibins - 1;
	    //  //if (iphi < 0 || iphi >= nphibins) std::cout << iphi << endl;
            //  dN2dphidy_ccbar_mixed_ALICE_v2[iphi] += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
	    //  ccbar_mixed_pairs += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
            //}
	  }
	}

	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy>=0 && iy<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdy_ccbar[iy] += (double)nybins/(2.*DNDY_MAX);
	    double initialmag = itsIP.px()*itsIP.px() + itsIP.py()*itsIP.py() + itsIP.pz()*itsIP.pz();
	    double eopdoti = itsIP.px()*itsEOPP.px() + itsIP.py()*itsEOPP.py() + itsIP.pz()*itsEOPP.pz();
	    double EOPdelpt2 = pow((eopdoti/initialmag)*itsIP.px() - itsEOPP.px(),2.);
	    EOPdelpt2 += pow((eopdoti/initialmag)*itsIP.py() - itsEOPP.py(),2.);
	    EOPdelpt2 += pow((eopdoti/initialmag)*itsIP.pz() - itsEOPP.pz(),2.);
            delpt2_preeq[iy] += EOPdelpt2;
	    double doti = itsIP.px()*itsP.px() + itsIP.py()*itsP.py() + itsIP.pz()*itsP.pz();
	    double Fdelpt2 = pow((doti/initialmag)*itsIP.px() - itsP.px(),2.);
	    Fdelpt2 += pow((doti/initialmag)*itsIP.py() - itsP.py(),2.);
	    Fdelpt2 += pow((doti/initialmag)*itsIP.pz() - itsP.pz(),2.);
            delpt2[iy] += Fdelpt2;
	    double y0 = 0.5*log((itsIP.e() + itsIP.pz())/(itsIP.e() - itsIP.pz()));
	    double yf = 0.5*log((itsP.e() + itsP.pz())/(itsP.e() - itsP.pz()));
	    dely2[iy] += (yf-y0)*(yf-y0);
	    //Now calculate this ignoring pz
	    initialmag = itsIP.px()*itsIP.px() + itsIP.py()*itsIP.py();
	    eopdoti = itsIP.px()*itsEOPP.px() + itsIP.py()*itsEOPP.py();
	    EOPdelpt2 = pow((eopdoti/initialmag)*itsIP.px() - itsEOPP.px(),2.);
	    EOPdelpt2 += pow((eopdoti/initialmag)*itsIP.py() - itsEOPP.py(),2.);
            delpt2_nopz_preeq[iy] += EOPdelpt2;
	    doti = itsIP.px()*itsP.px() + itsIP.py()*itsP.py();
	    Fdelpt2 = pow((doti/initialmag)*itsIP.px() - itsP.px(),2.);
	    Fdelpt2 += pow((doti/initialmag)*itsIP.py() - itsP.py(),2.);
            delpt2_nopz[iy] += Fdelpt2;
	  }
	}
	int ipseudorap = (int)( (pseudorap+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(ipseudorap>=0 && ipseudorap<nybins){
	  if(abs(plist[0]->at(ie).id())==4){
	    dNdpseudorap_ccbar[ipseudorap] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX_CMS && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ccbar_CMS[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ccbar_CMS[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_bbbar_CMS[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_bbbar_CMS[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }
	if(fabs(y)<Y_MAX_ALICE_v2 && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ccbar_ALICE_v2[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ccbar_ALICE_v2[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_bbbar_ALICE_v2[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_bbbar_ALICE_v2[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }
	if(fabs(y)<Y_MAX_ALICE_RAA && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(plist[0]->at(ie).id())==4){
	    dN2piptdptdy_ccbar_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(abs(plist[0]->at(ie).id())==5){
	    dN2piptdptdy_bbbar_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }
      }
	  
      eventSize = martini.pythia.event.size();

      for(int ie=0; ie<eventSize; ie++){
	Vec4 itsP = martini.pythia.event[ie].p();
	double pt = sqrt(itsP.px()*itsP.px() + itsP.py()*itsP.py());
	double itsM = itsP.mCalc();
	double pz = itsP.pz();
	//rapidity is more appropriate here:
	//double y = 0.5*log( (sqrt(pt*pt+pz*pz) + pz)/(sqrt(pt*pt+pz*pz) - pz) );
	double y = 0.5*log( (sqrt(pt*pt+pz*pz + itsM*itsM) + pz)/(sqrt(pt*pt+pz*pz + itsM*itsM) - pz) );
	double cos2Phi = 0.;
	double sin2Phi = 0.;
	if(pt > 1e-5){
	  cos2Phi = (itsP.px()*itsP.px() - itsP.py()*itsP.py())/(pt*pt);
	  sin2Phi = 2.*(itsP.px()*itsP.py())/(pt*pt);
	}

	if (martini.pythia.event[ie].id()==411)
        {
	  for (int je=0; je<eventSize; je++){
            if (martini.pythia.event[je].id() != -411) continue;
	    Vec4 itsPj = martini.pythia.event[je].p();
	    double ptj = sqrt(itsPj.px()*itsPj.px() + itsPj.py()*itsPj.py());
	    double itsMj = itsPj.mCalc();
            double pzj = itsPj.pz();
            //rapidity is more appropriate here:
            double yj         = 0.5*log( (sqrt(ptj*ptj+pzj*pzj+itsMj*itsMj) + pzj)/(sqrt(ptj*ptj+pzj*pzj+itsMj*itsMj) - pzj) );
	    if(fabs(y)<Y_MAX_dNdphi && fabs(yj)<Y_MAX_dNdphi){
	      double phi1 = atan(itsP.py()/itsP.px());
	      if (itsP.px() < 0) phi1 = phi1 + M_PI;
	      double phi2 = atan(itsPj.py()/itsPj.px());
	      if (itsPj.px() < 0) phi2 = phi2 + M_PI;
	      double delphi = phi1-phi2;
	      if (delphi < 0) delphi = - delphi;
	      if (delphi > M_PI) delphi = 2.*M_PI - delphi;
	      int iphi =(int)( delphi*(double)nphibins/(M_PI));
	      if (iphi == nphibins) iphi = nphibins - 1;
	      //if (iphi < 0 || iphi >= nphibins) std::cout << iphi << endl;
              dN2dphidy_DDbar[iphi] += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
	      DDbar_pairs += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
            }
	    //if(fabs(y)<Y_MAX_ALICE_v2 && yj<(Y_MAX_ALICE_v2+3.) && yj>(3.-Y_MAX_ALICE_v2)){
	    //  double phi1 = atan(itsP.py()/itsP.px());
	    //  if (itsP.px() < 0) phi1 = phi1 + M_PI;
	    //  double phi2 = atan(itsPj.py()/itsPj.px());
	    //  if (itsPj.px() < 0) phi2 = phi2 + M_PI;
	    //  double delphi = phi1-phi2;
	    //  if (delphi < 0) delphi = - delphi;
	    //  if (delphi > M_PI) delphi = 2.*M_PI - delphi;
	    //  int iphi =(int)( delphi*(double)nphibins/(M_PI));
	    //  if (iphi == nphibins) iphi = nphibins - 1;
	    //  //if (iphi < 0 || iphi >= nphibins) std::cout << iphi << endl;
            //  dN2dphidy_DDbar_mixed_ALICE_v2[iphi] += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
	    //  DDbar_mixed_pairs += 1.;//(double)nphibins/(M_PI*2.*Y_MAX_ALICE_v2);
            //}
	  }
	}
	//The dN/dy histograms:
	int iy = (int)( (y+DNDY_MAX)*(double)nybins/(2.*DNDY_MAX) );
	if(iy >= 0 && iy < nybins){
	  if(abs(martini.pythia.event[ie].id())==411){
	    dNdy_DDbar[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	  if(martini.pythia.event[ie].id()==443){
	    dNdy_JPsi[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	  if(martini.pythia.event[ie].id()==443000){
	    dNdy_ExcitedJPsi[iy] += (double)nybins/(2.*DNDY_MAX);
	  }
	}

	if(fabs(y)<Y_MAX_CMS && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(martini.pythia.event[ie].id())==411){
	    dN2piptdptdy_DDbar_CMS[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DDbar_CMS[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DDbar_CMS[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_CMS)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }

	if(fabs(y)<Y_MAX_ALICE_v2 && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(martini.pythia.event[ie].id())==411){
	    dN2piptdptdy_DDbar_ALICE_v2[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DDbar_ALICE_v2[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DDbar_ALICE_v2[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_v2)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }

	if(fabs(y)<Y_MAX_ALICE_RAA && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );
	  if(abs(martini.pythia.event[ie].id())==411){
	    dN2piptdptdy_DDbar_ALICE_RAA[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX_ALICE_RAA)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
        }

	if(fabs(y)<Y_MAX && pt<PT_MAX){
	  int ipt = (int)((pt/PT_MAX)*nptbins );

	  if(abs(martini.pythia.event[ie].id())==511){
	    dN2piptdptdy_DBDBbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_DBDBbar[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_DBDBbar[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	  }
	  if(martini.pythia.event[ie].id()==443){
	    dN2piptdptdy_JPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_JPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_JPsi[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    if(martini.pythia.event[ie].status()== 81){
	      dN2piptdptdy_OnlyDirectJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      sin2Phi_dN2piptdptdy_OnlyDirectJPsi[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    }
	    NJPsi += 1.;
	  }
	  if(martini.pythia.event[ie].id()==443000){
	    dN2piptdptdy_ExcitedJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    cos2Phi_dN2piptdptdy_ExcitedJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    sin2Phi_dN2piptdptdy_ExcitedJPsi[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NExcitedJPsi += 1.;
	    if(martini.pythia.event[ie].status()== 81){
	      dN2piptdptdy_OnlyDirectExcitedJPsi[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ipt] += cos2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	      sin2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ipt] += sin2Phi*(PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    }
	  }
	  if(martini.pythia.event[ie].id()==553){
	    dN2piptdptdy_Upsilon[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NUpsilon += 1.;
	  }
	  if(martini.pythia.event[ie].id()==553000){
	    dN2piptdptdy_ExcitedUpsilon[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NExcitedUpsilon += 1.;
	  }
	  if(martini.pythia.event[ie].id()==543){
	    dN2piptdptdy_bcbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NBc += 1.;
	  }
	  if(martini.pythia.event[ie].id()==543000){
	    dN2piptdptdy_Excitedbcbar[ipt] += (PT_MAX/(double)nptbins)*(0.5/Y_MAX)/(2.*M_PI*PT_MAX*((double)ipt+0.5)/(double)nptbins );
	    NBc += 1.;
	  }
	}
      }
    }

  //Output the properly normalized numbers of quarkonia:
  NJPsi /= (double)numJets;
  NExcitedJPsi /= (double)numJets;
  NUpsilon /= (double)numJets;
  NExcitedUpsilon /= (double)numJets;
  NBc /= (double)numJets;
  NExcitedBc /= (double)numJets;

  cout << "The average number of JPsi particles per collision = " << NJPsi << endl;
  cout << "The average number of excited charmonia per collision = " << NExcitedJPsi << endl;
  cout << "The average number of Upsilon particles per collision = " << NUpsilon << endl;
  cout << "The average number of excited bottomonia per collision = " << NExcitedUpsilon << endl;
  cout << "The average number of B_c mesons per collision = " << NBc << endl;
  cout << "The average number of excited B_c mesons per collision = " << NExcitedBc << endl;

  stringstream s;
  s << argv[1];
  
  if(outputting){

    string ccbarBeforeEvolutionnameCMS = directory_name+"dN2piptdpt_ccbar_beforeEvolution_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOutCMS;
    ccbarBeforeEvolutionOutCMS.open(ccbarBeforeEvolutionnameCMS.c_str());

    string bbbarBeforeEvolutionnameCMS = directory_name+"dN2piptdpt_bbbar_beforeEvolution_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarBeforeEvolutionOutCMS;
    bbbarBeforeEvolutionOutCMS.open(bbbarBeforeEvolutionnameCMS.c_str());

    string ccbarBeforeEvolutionnameALICEv2 = directory_name+"dN2piptdpt_ccbar_beforeEvolution_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOutALICEv2;
    ccbarBeforeEvolutionOutALICEv2.open(ccbarBeforeEvolutionnameALICEv2.c_str());

    string bbbarBeforeEvolutionnameALICEv2 = directory_name+"dN2piptdpt_bbbar_beforeEvolution_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarBeforeEvolutionOutALICEv2;
    bbbarBeforeEvolutionOutALICEv2.open(bbbarBeforeEvolutionnameALICEv2.c_str());

    string ccbarBeforeEvolutionnameALICERAA = directory_name+"dN2piptdpt_ccbar_beforeEvolution_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutionOutALICERAA;
    ccbarBeforeEvolutionOutALICERAA.open(ccbarBeforeEvolutionnameALICERAA.c_str());

    string bbbarBeforeEvolutionnameALICERAA = directory_name+"dN2piptdpt_bbbar_beforeEvolution_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarBeforeEvolutionOutALICERAA;
    bbbarBeforeEvolutionOutALICERAA.open(bbbarBeforeEvolutionnameALICERAA.c_str());

    string ccbardphiBeforeEvolutionname = directory_name+"dN2dphidy_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardphiBeforeEvolutionOut;
    ccbardphiBeforeEvolutionOut.open(ccbardphiBeforeEvolutionname.c_str());

//    string ccbarBeforeHydroname = directory_name+"dN2piptdpt_ccbar_beforeHydro_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarBeforeHydroOut;
//    ccbarBeforeHydroOut.open(ccbarBeforeHydroname.c_str());
//
//    string bbbarBeforeHydroname = directory_name+"dN2piptdpt_bbbar_beforeHydro_"+HQ_name+"_"+s.str()+".dat";
//    ofstream bbbarBeforeHydroOut;
//    bbbarBeforeHydroOut.open(bbbarBeforeHydroname.c_str());
//
//    string ccbarBeforeEvolutionLargeptname = directory_name+"dN2piptdpt_ccbar_beforeEvolution_largept_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarBeforeEvolutionLargeptOut;
//    ccbarBeforeEvolutionLargeptOut.open(ccbarBeforeEvolutionLargeptname.c_str());

//    string ccbarBeforeHydroLargeptname = directory_name+"dN2piptdpt_ccbar_beforeHydro_largept_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarBeforeHydroLargeptOut;
//    ccbarBeforeHydroLargeptOut.open(ccbarBeforeHydroLargeptname.c_str());
//
//    string ccbarLargeptname = directory_name+"dN2piptdpt_ccbar_largept_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarLargeptOut;
//    ccbarLargeptOut.open(ccbarLargeptname.c_str());

    string ccbarnameCMS = directory_name+"dN2piptdpt_ccbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarOutCMS;
    ccbarOutCMS.open(ccbarnameCMS.c_str());
    string cos2Phi_ccbarnameCMS = directory_name+"cos2Phi_dN2piptdpt_ccbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ccbarOutCMS;
    cos2Phi_ccbarOutCMS.open(cos2Phi_ccbarnameCMS.c_str());

    string bbbarnameCMS = directory_name+"dN2piptdpt_bbbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarOutCMS;
    bbbarOutCMS.open(bbbarnameCMS.c_str());
    string cos2Phi_bbbarnameCMS = directory_name+"cos2Phi_dN2piptdpt_bbbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_bbbarOutCMS;
    cos2Phi_bbbarOutCMS.open(cos2Phi_bbbarnameCMS.c_str());

    string DDbarnameCMS = directory_name+"dN2piptdpt_DDbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbarOutCMS;
    DDbarOutCMS.open(DDbarnameCMS.c_str());
    string cos2Phi_DDbarnameCMS = directory_name+"cos2Phi_dN2piptdpt_DDbar_CMS_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DDbarOutCMS;
    cos2Phi_DDbarOutCMS.open(cos2Phi_DDbarnameCMS.c_str());

    string ccbarnameALICEv2 = directory_name+"dN2piptdpt_ccbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarOutALICEv2;
    ccbarOutALICEv2.open(ccbarnameALICEv2.c_str());
    string cos2Phi_ccbarnameALICEv2 = directory_name+"cos2Phi_dN2piptdpt_ccbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ccbarOutALICEv2;
    cos2Phi_ccbarOutALICEv2.open(cos2Phi_ccbarnameALICEv2.c_str());

    string ccbardphiname = directory_name+"dN2dphidy_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardphiOut;
    ccbardphiOut.open(ccbardphiname.c_str());

    string bbbarnameALICEv2 = directory_name+"dN2piptdpt_bbbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarOutALICEv2;
    bbbarOutALICEv2.open(bbbarnameALICEv2.c_str());
    string cos2Phi_bbbarnameALICEv2 = directory_name+"cos2Phi_dN2piptdpt_bbbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_bbbarOutALICEv2;
    cos2Phi_bbbarOutALICEv2.open(cos2Phi_bbbarnameALICEv2.c_str());

    string DDbarnameALICEv2 = directory_name+"dN2piptdpt_DDbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbarOutALICEv2;
    DDbarOutALICEv2.open(DDbarnameALICEv2.c_str());
    string cos2Phi_DDbarnameALICEv2 = directory_name+"cos2Phi_dN2piptdpt_DDbar_ALICEv2_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DDbarOutALICEv2;
    cos2Phi_DDbarOutALICEv2.open(cos2Phi_DDbarnameALICEv2.c_str());

    string DDbardphiname = directory_name+"dN2dphidy_DDbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbardphiOut;
    DDbardphiOut.open(DDbardphiname.c_str());

    string ccbarnameALICERAA = directory_name+"dN2piptdpt_ccbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarOutALICERAA;
    ccbarOutALICERAA.open(ccbarnameALICERAA.c_str());
    string cos2Phi_ccbarnameALICERAA = directory_name+"cos2Phi_dN2piptdpt_ccbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ccbarOutALICERAA;
    cos2Phi_ccbarOutALICERAA.open(cos2Phi_ccbarnameALICERAA.c_str());

    string bbbarnameALICERAA = directory_name+"dN2piptdpt_bbbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream bbbarOutALICERAA;
    bbbarOutALICERAA.open(bbbarnameALICERAA.c_str());
    string cos2Phi_bbbarnameALICERAA = directory_name+"cos2Phi_dN2piptdpt_bbbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_bbbarOutALICERAA;
    cos2Phi_bbbarOutALICERAA.open(cos2Phi_bbbarnameALICERAA.c_str());

    string DDbarnameALICERAA = directory_name+"dN2piptdpt_DDbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbarOutALICERAA;
    DDbarOutALICERAA.open(DDbarnameALICERAA.c_str());
    string cos2Phi_DDbarnameALICERAA = directory_name+"cos2Phi_dN2piptdpt_DDbar_ALICERAA_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DDbarOutALICERAA;
    cos2Phi_DDbarOutALICERAA.open(cos2Phi_DDbarnameALICERAA.c_str());

    string DBDBbarname = directory_name+"dN2piptdpt_DBDBbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream DBDBbarOut;
    DBDBbarOut.open(DBDBbarname.c_str());
    string cos2Phi_DBDBbarname = directory_name+"cos2Phi_dN2piptdpt_DBDBbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_DBDBbarOut;
    cos2Phi_DBDBbarOut.open(cos2Phi_DBDBbarname.c_str());

    string JPsiname = directory_name+"dN2piptdpt_JPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream JPsiOut;
    JPsiOut.open(JPsiname.c_str());
    string cos2Phi_JPsiname = directory_name+"cos2Phi_dN2piptdpt_JPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_JPsiOut;
    cos2Phi_JPsiOut.open(cos2Phi_JPsiname.c_str());

    string OnlyDirectJPsiname = directory_name+"dN2piptdpt_OnlyDirectJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream OnlyDirectJPsiOut;
    OnlyDirectJPsiOut.open(OnlyDirectJPsiname.c_str());
    string cos2Phi_OnlyDirectJPsiname = directory_name+"cos2Phi_dN2piptdpt_OnlyDirectJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_OnlyDirectJPsiOut;
    cos2Phi_OnlyDirectJPsiOut.open(cos2Phi_OnlyDirectJPsiname.c_str());

    string ExcitedJPsiname = directory_name+"dN2piptdpt_ExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedJPsiOut;
    ExcitedJPsiOut.open(ExcitedJPsiname.c_str());
    string cos2Phi_ExcitedJPsiname = directory_name+"cos2Phi_dN2piptdpt_ExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_ExcitedJPsiOut;
    cos2Phi_ExcitedJPsiOut.open(cos2Phi_ExcitedJPsiname.c_str());

    string OnlyDirectExcitedJPsiname = directory_name+"dN2piptdpt_OnlyDirectExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream OnlyDirectExcitedJPsiOut;
    OnlyDirectExcitedJPsiOut.open(OnlyDirectExcitedJPsiname.c_str());
    string cos2Phi_OnlyDirectExcitedJPsiname = directory_name+"cos2Phi_dN2piptdpt_OnlyDirectExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream cos2Phi_OnlyDirectExcitedJPsiOut;
    cos2Phi_OnlyDirectExcitedJPsiOut.open(cos2Phi_OnlyDirectExcitedJPsiname.c_str());

    string Upsilonname = directory_name+"dN2piptdpt_Upsilon_"+HQ_name+"_"+s.str()+".dat";
    ofstream UpsilonOut;
    UpsilonOut.open(Upsilonname.c_str());

    string ExcitedUpsilonname = directory_name+"dN2piptdpt_ExcitedUpsilon_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedUpsilonOut;
    ExcitedUpsilonOut.open(ExcitedUpsilonname.c_str());

    string bcbarname = directory_name+"dN2piptdpt_bcbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream bcbarOut;
    bcbarOut.open(bcbarname.c_str());

    string Excitedbcbarname = directory_name+"dN2piptdpt_Excitedbcbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedbcbarOut;
    ExcitedbcbarOut.open(Excitedbcbarname.c_str());
    
    string PrehydroTname = directory_name+"PrehydroT_"+HQ_name+"_"+s.str()+".dat";
    ofstream PrehydroTOut;
    PrehydroTOut.open(PrehydroTname.c_str());
    
    string hydroTname = directory_name+"hydroT_"+HQ_name+"_"+s.str()+".dat";
    ofstream hydroTOut;
    hydroTOut.open(hydroTname.c_str());
    
    //for(int ip=0; ip<largeptbins; ip++){
    //  ccbarBeforeEvolutionLargeptOut << midpt[ip] << " " << dN2piptdptdy_ccbar_beforeEvolution_largept[ip]/(2.*event_counter) << endl;
    //  ccbarBeforeHydroLargeptOut << midpt[ip] << " " << dN2piptdptdy_ccbar_beforeHydro_largept[ip]/(2.*event_counter) << endl;
    //  ccbarLargeptOut << midpt[ip] << " " << dN2piptdptdy_ccbar_largept[ip]/(2.*event_counter) << endl;
    //}

    for(int iphi=0; iphi<nphibins; iphi++){
      ccbardphiBeforeEvolutionOut << M_PI*((double)iphi+0.5)/(double)nphibins << " " << dN2dphidy_ccbar_beforeEvolution[iphi]/ccbar_pairs_beforeevolution << endl;
      ccbardphiOut                << M_PI*((double)iphi+0.5)/(double)nphibins << " " << dN2dphidy_ccbar[iphi]/ccbar_pairs                                 << endl;
      DDbardphiOut                << M_PI*((double)iphi+0.5)/(double)nphibins << " " << dN2dphidy_DDbar[iphi]/DDbar_pairs                                 << endl;
      //cout << ccbar_pairs_beforeevolution << " " << ccbar_mixed_pairs_beforeevolution << " " << ccbar_pairs << " " << ccbar_mixed_pairs << " " << DDbar_pairs << " " << DDbar_mixed_pairs << endl;
    }

    for(int ip=0; ip<nptbins; ip++){
      ccbarBeforeEvolutionOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeEvolution_CMS[ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeEvolution_CMS[ip]/(2.*event_counter) << endl;

      ccbarBeforeEvolutionOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeEvolution_ALICE_v2[ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeEvolution_ALICE_v2[ip]/(2.*event_counter) << endl;

      ccbarBeforeEvolutionOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeEvolution_ALICE_RAA[ip]/(2.*event_counter) << endl;
      bbbarBeforeEvolutionOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeEvolution_ALICE_RAA[ip]/(2.*event_counter) << endl;

//      ccbarBeforeHydroOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_beforeHydro[ip]/(2.*event_counter) << endl;
//      bbbarBeforeHydroOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_beforeHydro[ip]/(2.*event_counter) << endl;

      ccbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_CMS[ip]/(2.*event_counter) << endl;
      bbbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_CMS[ip]/(2.*event_counter) << endl;
      cos2Phi_ccbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ccbar_CMS[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_ccbar_CMS[ip]/(2.*event_counter) << endl;
      cos2Phi_bbbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_bbbar_CMS[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_bbbar_CMS[ip]/(2.*event_counter) << endl;

      ccbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_ALICE_v2[ip]/(2.*event_counter) << endl;
      bbbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_ALICE_v2[ip]/(2.*event_counter) << endl;
      cos2Phi_ccbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ccbar_ALICE_v2[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_ccbar_ALICE_v2[ip]/(2.*event_counter) << endl;
      cos2Phi_bbbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_bbbar_ALICE_v2[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_bbbar_ALICE_v2[ip]/(2.*event_counter) << endl;

      ccbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ccbar_ALICE_RAA[ip]/(2.*event_counter) << endl;
      bbbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bbbar_ALICE_RAA[ip]/(2.*event_counter) << endl;
      cos2Phi_ccbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_ccbar_ALICE_RAA[ip]/(2.*event_counter) << endl;
      cos2Phi_bbbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_bbbar_ALICE_RAA[ip]/(2.*event_counter) << endl;

      DDbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar_CMS[ip]/(2.*event_counter) << endl;
      cos2Phi_DDbarOutCMS << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DDbar_CMS[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DDbar_CMS[ip]/(2.*event_counter) << endl;

      DDbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar_ALICE_v2[ip]/(2.*event_counter) << endl;
      cos2Phi_DDbarOutALICEv2 << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DDbar_ALICE_v2[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DDbar_ALICE_v2[ip]/(2.*event_counter) << endl;

      DDbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DDbar_ALICE_RAA[ip]/(2.*event_counter) << endl;
      cos2Phi_DDbarOutALICERAA << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DDbar_ALICE_RAA[ip]/(2.*event_counter) << endl;

      DBDBbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << endl;
      cos2Phi_DBDBbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << " " << sin2Phi_dN2piptdptdy_DBDBbar[ip]/(2.*event_counter) << endl;

      JPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_JPsi[ip]/(double)event_counter << endl;
      ExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << endl;
      cos2Phi_JPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_JPsi[ip]/(double)event_counter << " " << sin2Phi_dN2piptdptdy_JPsi[ip]/(double)event_counter << endl;
      cos2Phi_ExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << " " << sin2Phi_dN2piptdptdy_ExcitedJPsi[ip]/(double)event_counter << endl;
      OnlyDirectJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_OnlyDirectJPsi[ip]/(double)event_counter << endl;
      OnlyDirectExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_OnlyDirectExcitedJPsi[ip]/(double)event_counter << endl;
      cos2Phi_OnlyDirectJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_OnlyDirectJPsi[ip]/(double)event_counter << " " << sin2Phi_dN2piptdptdy_OnlyDirectJPsi[ip]/(double)event_counter << endl;
      cos2Phi_OnlyDirectExcitedJPsiOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << cos2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip]/(double)event_counter << " " << sin2Phi_dN2piptdptdy_OnlyDirectExcitedJPsi[ip]/(double)event_counter << endl;

      UpsilonOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_Upsilon[ip]/(double)event_counter << endl;
      ExcitedUpsilonOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_ExcitedUpsilon[ip]/(double)event_counter << endl;

      bcbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_bcbar[ip]/(double)event_counter << endl;
      ExcitedbcbarOut << PT_MAX*((double)ip+0.5)/(double)nptbins << " " << dN2piptdptdy_Excitedbcbar[ip]/(double)event_counter << endl;
    }
  
    for (int ind = 0; ind < 100; ind++){
      //if (prehydro_counter[ind] == 0) break;
      PrehydroTOut << ind << " " << prehydro_T[ind] << " " << prehydro_counter[ind] << endl;
    }

    for (int ind = 0; ind < 1000; ind++){
      if (hydro_counter[ind] == 0) break;
      hydroTOut << ind << " " << hydro_T[ind] << " " << hydro_counter[ind] << endl;
    }

    string ccbarBeforeEvolutiondNdyname = directory_name+"dNdy_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdyOut;
    ccbarBeforeEvolutiondNdyOut.open(ccbarBeforeEvolutiondNdyname.c_str());

    string ccbarBeforeEvolutiondNdpseudorapname = directory_name+"dNdpseudorap_ccbar_beforeEvolution_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbarBeforeEvolutiondNdpseudorapOut;
    ccbarBeforeEvolutiondNdpseudorapOut.open(ccbarBeforeEvolutiondNdpseudorapname.c_str());

//    string ccbarBeforeHydrodNdyname = directory_name+"dNdy_ccbar_beforeHydro_"+HQ_name+"_"+s.str()+".dat";
//    ofstream ccbarBeforeHydrodNdyOut;
//    ccbarBeforeHydrodNdyOut.open(ccbarBeforeHydrodNdyname.c_str());

    string ccbardNdyname = directory_name+"dNdy_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdyOut;
    ccbardNdyOut.open(ccbardNdyname.c_str());

    string ccbardNdpseudorapname = directory_name+"dNdpseudorap_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardNdpseudorapOut;
    ccbardNdpseudorapOut.open(ccbardNdpseudorapname.c_str());

    string ccbardelpt2preeqname = directory_name+"delpT2dy_ccbar_preeq_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardelpt2preeqOut;
    ccbardelpt2preeqOut.open(ccbardelpt2preeqname.c_str());

    string ccbardelpt2name = directory_name+"delpT2dy_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardelpt2Out;
    ccbardelpt2Out.open(ccbardelpt2name.c_str());

    string ccbardely2name = directory_name+"dely2dy_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardely2Out;
    ccbardely2Out.open(ccbardely2name.c_str());

    string ccbardelpt2nopzpreeqname = directory_name+"delpT2dy_nopz_ccbar_preeq_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardelpt2nopzpreeqOut;
    ccbardelpt2nopzpreeqOut.open(ccbardelpt2nopzpreeqname.c_str());

    string ccbardelpt2nopzname = directory_name+"delpT2dy_nopz_ccbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream ccbardelpt2nopzOut;
    ccbardelpt2nopzOut.open(ccbardelpt2nopzname.c_str());

    string DDbardNdyname = directory_name+"dNdy_DDbar_"+HQ_name+"_"+s.str()+".dat";
    ofstream DDbardNdyOut;
    DDbardNdyOut.open(DDbardNdyname.c_str());

    string JPsidNdyname = directory_name+"dNdy_JPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream JPsidNdyOut;
    JPsidNdyOut.open(JPsidNdyname.c_str());

    string ExcitedJPsidNdyname = directory_name+"dNdy_ExcitedJPsi_"+HQ_name+"_"+s.str()+".dat";
    ofstream ExcitedJPsidNdyOut;
    ExcitedJPsidNdyOut.open(ExcitedJPsidNdyname.c_str());
    
    for(int iy=0; iy<nybins; iy++){
      ccbarBeforeEvolutiondNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar_beforeEvolution[iy]/(2.*(double)event_counter) << endl;
      ccbarBeforeEvolutiondNdpseudorapOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdpseudorap_ccbar_beforeEvolution[iy]/(2.*(double)event_counter) << endl;
    //  ccbarBeforeHydrodNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar_beforeHydro[iy]/(2.*(double)event_counter) << endl;
      ccbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ccbar[iy]/(2.*(double)event_counter) << endl;
      ccbardNdpseudorapOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdpseudorap_ccbar[iy]/(2.*(double)event_counter) << endl;
      ccbardelpt2preeqOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << (delpt2_preeq[iy]/dNdy_ccbar[iy])*(double)nybins/(2.*DNDY_MAX) << endl;
      ccbardelpt2Out << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << (delpt2[iy]/dNdy_ccbar[iy])*(double)nybins/(2.*DNDY_MAX) << endl;
      ccbardely2Out << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << (dely2[iy]/dNdy_ccbar[iy])*(double)nybins/(2.*DNDY_MAX) << endl;
      ccbardelpt2nopzpreeqOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << (delpt2_nopz_preeq[iy]/dNdy_ccbar[iy])*(double)nybins/(2.*DNDY_MAX) << endl;
      ccbardelpt2nopzOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << (delpt2_nopz[iy]/dNdy_ccbar[iy])*(double)nybins/(2.*DNDY_MAX) << endl;
      DDbardNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_DDbar[iy]/(2.*(double)event_counter) << endl;
      JPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_JPsi[iy]/(2.*(double)event_counter) << endl;
      ExcitedJPsidNdyOut << (2.*((double)iy+0.5)/(double)nybins-1.)*DNDY_MAX << " " << dNdy_ExcitedJPsi[iy]/(2.*(double)event_counter) << endl;
    }

    PrehydroTOut.close();    
    hydroTOut.close();   

    ccbardelpt2preeqOut.close();
    ccbardelpt2Out.close();
    ccbardely2Out.close();
    ccbardelpt2nopzpreeqOut.close();
    ccbardelpt2nopzOut.close();

    ccbarBeforeEvolutionOutCMS.close();
    bbbarBeforeEvolutionOutCMS.close();

    ccbarBeforeEvolutionOutALICEv2.close();
    bbbarBeforeEvolutionOutALICEv2.close();

    ccbarBeforeEvolutionOutALICERAA.close();
    bbbarBeforeEvolutionOutALICERAA.close();

    ccbardphiBeforeEvolutionOut.close();
    ccbardphiOut.close();
    DDbardphiOut.close();
//    ccbarBeforeHydroOut.close();

//    bbbarBeforeHydroOut.close();

//    ccbarBeforeEvolutionLargeptOut.close();
//    ccbarBeforeHydroLargeptOut.close();
//    ccbarLargeptOut.close();

    ccbarOutCMS.close();
    bbbarOutCMS.close();

    ccbarOutALICEv2.close();
    bbbarOutALICEv2.close();

    ccbarOutALICERAA.close();
    bbbarOutALICERAA.close();

    cos2Phi_ccbarOutCMS.close();
    cos2Phi_bbbarOutCMS.close();

    cos2Phi_ccbarOutALICEv2.close();
    cos2Phi_bbbarOutALICEv2.close();

    cos2Phi_ccbarOutALICERAA.close();
    cos2Phi_bbbarOutALICERAA.close();

    DDbarOutCMS.close();
    cos2Phi_DDbarOutCMS.close();

    DDbarOutALICEv2.close();
    cos2Phi_DDbarOutALICEv2.close();

    DDbarOutALICERAA.close();
    cos2Phi_DDbarOutALICERAA.close();

    DBDBbarOut.close();
    cos2Phi_DBDBbarOut.close();

    JPsiOut.close();
    cos2Phi_JPsiOut.close();
    OnlyDirectJPsiOut.close();
    cos2Phi_OnlyDirectJPsiOut.close();

    ExcitedJPsiOut.close();
    cos2Phi_ExcitedJPsiOut.close();
    OnlyDirectExcitedJPsiOut.close();
    cos2Phi_OnlyDirectExcitedJPsiOut.close();

    UpsilonOut.close();

    ExcitedUpsilonOut.close();

    bcbarOut.close();

    ExcitedbcbarOut.close();
    
    ccbarBeforeEvolutiondNdyOut.close();
    ccbarBeforeEvolutiondNdpseudorapOut.close();
    //ccbarBeforeHydrodNdyOut.close();
    ccbardNdyOut.close();
    ccbardNdpseudorapOut.close();
    DDbardNdyOut.close();
    JPsidNdyOut.close();
    ExcitedJPsidNdyOut.close();

  }

  delete plist[0];
  delete plist[1];
  delete plist;  
  delete plistInitial;
}
