#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <sstream>

using namespace std;

inline string stringify(int x);

double edt,elt,clt,ps,hms_rate,hms34,hms_tof,hms_prlo,hms_trackeff,bcm1,bcm2,q1,q2,eb,p,theta;
int target,rstart,rend,pol,com1,com2,events;
int i,fcheck;

int main () {
  string t1,t2,t3,t4,t5,t6,t7,t8,syspar,runnum,scname,epicname;
  string scpath="/w/work5501/inclusive/replay/pass2/scalers/";
  ifstream scfile;
  fstream sctmpfile;
  ofstream outfile;
  
  outfile.open("rosen07_runlist.dat");

  cout << "runstart:";
  cin >> rstart;
  cout << "run end:";
  cin >> rend;  
  cout <<setw(6)<<right<<"!run#"
       <<setw(7)<<"Eb"
       <<setw(4)<<"Tar"
       <<setw(8)<<"theta"
       <<setw(8)<<"Ep"
       <<setw(4)<<"Pol"
       <<setw(10)<<"Events"
       <<setw(5)<<"com1"
       <<setw(5)<<"com2"
       <<setw(7)<<"ps"
       <<setw(12)<<"bcm1"
       <<setw(12)<<"bcm2"
       <<setw(12)<<"q1"
       <<setw(12)<<"q2"
       <<setw(7)<<"clt"
       <<setw(7)<<"elt"
       <<setw(8)<<"h_treff"
       <<setw(8)<<"h_tof"
       <<setw(8)<<"hms34"
       <<setw(8)<<"h_prlo"
       <<setw(11)<<"h_rate"<<endl;
  cout <<setw(160)<<setfill('-')<<"!"<<endl;
  outfile <<setw(6)<<right<<"run#"
	  <<setw(7)<<"Eb"
	  <<setw(4)<<"Tar"
	  <<setw(8)<<"theta"
	  <<setw(8)<<"Ep"
	  <<setw(4)<<"Pol"
	  <<setw(10)<<"Events"
	  <<setw(5)<<"com1"
	  <<setw(5)<<"com2"
	  <<setw(7)<<"ps"
	  <<setw(12)<<"bcm1"
	  <<setw(12)<<"bcm2"
	  <<setw(12)<<"q1"
	  <<setw(12)<<"q2"
	  <<setw(7)<<"clt"
	  <<setw(7)<<"elt"
	  <<setw(8)<<"h_treff"
	  <<setw(8)<<"h_tof"
	  <<setw(8)<<"hms34"
	  <<setw(8)<<"h_prlo"
	  <<setw(11)<<"h_rate"<<endl;
  outfile <<setw(160)<<setfill('-')<<"!"<<endl;


  for(i=rstart;i<=rend;i++){
 
    runnum = stringify(i);
    scname = scpath+"gen"+runnum+".txt";
    epicname = scpath+"epics"+runnum+".txt";
    scfile.open(scname.c_str());
    if (scfile.is_open()) {
      //cout << scname << endl;
      /********* Get Scaler Info*******************/
      syspar = "./get_scaler_info.sh "+runnum+" > scals.tmp";
      system("rm scals.tmp");
      system(syspar.c_str());
      sctmpfile.open ("scals.tmp");
      sctmpfile >> target>>eb>>p>>theta>>ps>>bcm1>>bcm2>>q1>>q2>>clt>>elt>>hms_trackeff>>hms_tof>>hms34>>hms_prlo>>hms_rate;
      sctmpfile.close();

      /********* Get number of events(Hcleantrack)**************/
      syspar = "grep 'sync hELCLEAN =' "+scname+" > scals.tmp";
      system("rm scals.tmp");
      system(syspar.c_str());
      sctmpfile.open ("scals.tmp");
      sctmpfile >> t1 >> t2 >> t3 >> events;
      sctmpfile.close();
      /*********************************************************/

      /********* Get HMS Polarity           ********************/
      syspar = "grep 'ecDI_POLARITY' "+epicname+" > scals.tmp";
      system("rm scals.tmp");
      system(syspar.c_str());
      sctmpfile.open ("scals.tmp");
      sctmpfile >> t1 >> pol;
      sctmpfile.close();
      /*********************************************************/

      cout.fill (' ');
      cout.width (10);

      cout <<setw(6)<<right<<runnum
	   <<setw(7)<<fixed<<setprecision(3)<<eb
	   <<setw(4)<<fixed<<setprecision(0)<<target
	   <<setw(8)<<fixed<<setprecision(2)<<theta
	   <<setw(8)<<fixed<<setprecision(4)<<p
	   <<setw(4)<<fixed<<setprecision(0)<<pol
	   <<setw(10)<<fixed<<setprecision(0)<<events
	   <<setw(5)<<fixed<<setprecision(0)<<com1
	   <<setw(5)<<fixed<<setprecision(0)<<com2
	   <<setw(7)<<fixed<<setprecision(3)<<ps
	   <<setw(12)<<fixed<<setprecision(2)<<bcm1
	   <<setw(12)<<fixed<<setprecision(2)<<bcm2
	   <<setw(12)<<fixed<<setprecision(2)<<q1
	   <<setw(12)<<fixed<<setprecision(2)<<q2
	   <<setw(7)<<fixed<<setprecision(4)<<clt
	   <<setw(7)<<elt
	   <<setw(8)<<hms_trackeff
	   <<setw(8)<<hms_tof
	   <<setw(8)<<hms34
	   <<setw(8)<<hms_prlo
	   <<setw(11)<<hms_rate<<endl;

      outfile.fill(' ');
      outfile.width(10);
      outfile <<setw(6)<<right<<runnum
	      <<setw(7)<<fixed<<setprecision(3)<<eb
	      <<setw(4)<<fixed<<setprecision(0)<<target
	      <<setw(8)<<fixed<<setprecision(2)<<theta
	      <<setw(8)<<fixed<<setprecision(4)<<p
	      <<setw(4)<<fixed<<setprecision(0)<<pol
	      <<setw(10)<<fixed<<setprecision(0)<<events
	      <<setw(5)<<fixed<<setprecision(0)<<com1
	      <<setw(5)<<fixed<<setprecision(0)<<com2
	      <<setw(7)<<fixed<<setprecision(3)<<ps
	      <<setw(12)<<fixed<<setprecision(2)<<bcm1
	      <<setw(12)<<fixed<<setprecision(2)<<bcm2
	      <<setw(12)<<fixed<<setprecision(2)<<q1
	      <<setw(12)<<fixed<<setprecision(2)<<q2
	      <<setw(7)<<fixed<<setprecision(4)<<clt
	      <<setw(7)<<elt
	      <<setw(8)<<hms_trackeff
	      <<setw(8)<<hms_tof
	      <<setw(8)<<hms34
	      <<setw(8)<<hms_prlo
	      <<setw(11)<<hms_rate<<endl;

      scfile.close();
    }   
  }
  outfile.close();

  return 0;
}

 inline string stringify(int x)
 {
   std::ostringstream o;
   o << x;
   return o.str();
 } 
