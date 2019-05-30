#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

float edt,elt,cdt,clt,ps,hms_rate,hms34,hms_tof,trkeff,hms_prlo,bcm1_I,bcm2_I,q1,q2,eb,p,theta;
int target;
int main () {
  string t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,syspar,runnum,scname,list;
  string scpath="/w/hallc/sane/anusha/replay/scalers/";
  fstream myfile;
  cin >> runnum;
  //runnum = "63929";
  scname=scpath+"gen"+runnum+".txt";
  //scname2=scpath+"gen"+runnum+".txt";

  /********* Get Electronic life time *******************/
  syspar = "grep 'Electronic  L.T.   :' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> elt >> t5;
  myfile.close();
  /**** elt = 1-edt/100;
  //cout << "E.Life Time:" << elt << endl;
  /*********************************************************/

  /********* Get Computer life time ************************/
  syspar = "grep 'Computer    D.T.   :' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> cdt;
  myfile.close();
  clt=1-cdt;
  //cout << "C.Life Time:" << clt << "C. Dead Time:"<< cdt << endl;
  /*********************************************************/

 /********* Get Target ************************************/
  syspar = "grep 'targ num' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> target >> t5;
  myfile.close();
  //  cout << "Target:" << target << endl;
  /*********************************************************/

  /********* Get Prescale factor  ************************************/
  syspar = "grep 'PS1 (HMS)   :' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> ps;
  myfile.close();
  //cout << "PS1(hms):" << ps << endl;
  /*********************************************************/

  /********* Get HMS rate  ************************************/
  //syspar = "grep 'Number of HMS      :' "+scname+" > scalers.tmp";
  syspar = "grep 'hSCIN          =' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> t4 >> hms_rate;
  hms_rate = hms_rate/1000; 
  myfile.close();
  // cout << "HMS rate:" << hms_rate << endl;
  /*********************************************************/

  /********* Get HMS 3/4 efficiency  ************************************/
  syspar = "grep 'HMS 3/4 effic  =' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> t4 >> hms34;
  myfile.close();
  //cout << "HMS 3/4 effic:" << hms34 << endl;
  /*********************************************************/

  /********* Get HMS TOF effic  ************************************/
  syspar = "grep 'HMS TOF effic  =' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> t4 >> hms_tof;
  myfile.close();
  //cout << "HMS tof effic:" << hms_tof << endl;
  /*********************************************************/

  /********* Get Tracking Eff  ************************************/
  syspar = "grep 'HMS efid effic  =' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> t4 >> trkeff;
  myfile.close();
  //cout << "HMS tof effic:" << hms_tof << endl;
  /*********************************************************/


  /********* Get HMS prlo  ************************************/
  syspar = "grep 'HMS prlo       =' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> hms_prlo;
  myfile.close();
  //  cout << "HMS prlo:" << hms_prlo << endl;
  /*********************************************************/

 /********* Get beam current & charge  ************************************/
  syspar = "grep 'Beamon BCM1' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> t4 >> t5 >> bcm1_I >> q1;
  myfile.close();
  syspar = "grep 'Beamon BCM2' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> t4 >> t5 >> bcm2_I >> q2;
  myfile.close();
  //cout << "bcm1:" << bcm1 << "q1:" << q1 << "bcm2:" << bcm2 << "q2:" << q2 << endl;
  /*********************************************************/

  /********* Get beam energy  ************************************/
  syspar = "grep 'E_beam' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> eb;
  myfile.close();
  //cout << "beam energy:" << eb << endl;
  /*********************************************************/

  /********* Get central momentum  ************************************/
  syspar = "grep 'P HMS' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> p;
  myfile.close();
  //cout << "central momentum:" << p << endl;
  /*********************************************************/

  /********* Get central angle  ************************************/
  syspar = "grep 'Theta HMS' "+scname+" > scalers.tmp";
  system("rm scalers.tmp");
  system(syspar.c_str());
  myfile.open ("scalers.tmp");
  myfile >> t1 >> t2 >> t3 >> theta;
  myfile.close();
  //cout << "central angle:" << theta << endl;
  /*********************************************************/

 
  cout << target << " " << eb << " " << p << " " << theta << " " << ps << " " << bcm1_I << " " << bcm2_I << " "; 
  cout << q1 << " " << q2 << " " << clt << " " << elt << " " << trkeff << " " << hms_tof << " " << hms34 << " " << hms_prlo << " " << hms_rate << endl;
 

 
  return 0;
}

