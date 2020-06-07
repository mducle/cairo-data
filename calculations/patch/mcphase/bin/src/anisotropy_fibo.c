/**************************************************************************
 *
 * anisotropy - program to calculate magnetic anisotropy
 *
 **************************************************************************/

#include<mcphas.h>

int verbose=1;
const char * filemode="w";
// for statistics 
int nofmaxloopDIV=0,nofmaxspinchangeDIV=0;
int successrate=0;
int nofcalls=0;
static const double PI_PHIINV = PI * 2. / (1. + sqrt(5.));
static const double TWOPI = PI * 2.;

#include "myev.h"
#include "mcphas_htcalc.c"
#include "mcphas_fecalc.c"
#include "mcphas_physpropcalc.c"

// main program
int main (int argc, char **argv)
{  std::clock_t startcputime = std::clock();
  char sipffilename[MAXNOFCHARINLINE];
  int im,l,nofsteps;
  int do_sipffile=0;
  int nofthreads=1;
  double z,u;
  double T,H;
  Vector xv(1,3);
  Vector yv(1,3);
  Vector h(1,3);
  
fprintf(stderr,"**************************************************************************\n");
fprintf(stderr,"*\n");
fprintf(stderr,"* anisotropy - program to calculate magnetic anisotropy \n");
fprintf(stderr,"*\n");
fprintf(stderr,"**************************************************************************\n\n");
  
// check command line
//T H xn yn zn nofsteps 
  if(argc<6){fprintf(stderr,"ERROR anisotropy: too few parameters\n");exit(EXIT_FAILURE);}
  T=strtod (argv[1], NULL); 
  H=strtod (argv[2], NULL); 
//Vector direction(1,3);
  //direction(1)=strtod (argv[3], NULL); 
  //direction(2)=strtod (argv[4], NULL); 
  //direction(3)=strtod (argv[5], NULL); 
  nofsteps=(int)strtod (argv[3], NULL); 
  
if (argc>=5)if (strcmp(argv[4],"-r")==0)
                 {do_sipffile=1;strcpy(sipffilename,argv[5]);
                 }

/*
 direction/=Norm(direction); // normalize
// now get r1 and r2 which are the basis vectors of the plane of rotation
 Vector x(1,3),r1(1,3),r2(1,3);x=0;x(1)=1;
 if (fabs(direction*x)>0.95){x=0;x(2)=1;}
 r1=x - (direction*x)*direction;
 r1/=Norm(r1);xproduct(r2,direction,r1);
 r2/=Norm(r2);
*/
  
 FILE * fout;fout=fopen_errchk("./results/anisotropy.out","w");
fprintf(fout,
"# output file of program: anisotropy @command\n"
"#! displayytext=M||H(mub)\n"
"#! displaytitle= Anisotropy plot\n"
"# phi(deg) theta(deg) T[K] |H|[T] Hx[T] Hy[T] Hz[T] |M|[mb] Mx[mb] My[mb] Mz[mb] MparallelH[mb]\n");

if(do_sipffile){
 jjjpar jjj(0,0,0,sipffilename,argc-9);jjj.save_sipf("./results/_");
 int nofcomponents=argc-9;
 Vector Hxc(1,nofcomponents);
 for(int j=1;j<=nofcomponents;++j)Hxc=(int)strtod (argv[j+8], NULL); 
 Hxc=0;h=0;Vector m(1,3);
 jjj.Icalc_parameter_storage_init(Hxc,h,T);
 // loop different H 
/*
 for(double az=0;az<2*PI-0.00001;az+=2*PI/nofsteps)
 {h=H*(cos(az)*r1+sin(az)*r2);
  double phi=PI/2;if(h(2)<0){phi=3*PI/2;}
  if(h(1)>0.001&&h(2)>=0){phi=atan(h(2)/h(1));}
  if(h(1)>0.001&&h(2)<0){phi=atan(h(2)/h(1))+2*PI;}
  if(h(1)<-0.001){phi=atan(h(2)/h(1))+PI;}
  double theta=acos(h(3)/H);*/
 for(int j=0; j<nofsteps; j++) {
  h(3) = 1. - (2. * j / nofsteps);
  double theta = acos(h(3));
  double phi = fmod(2. * j * PI_PHIINV, TWOPI);
  double stheta = sin(theta);
  h(1) = cos(phi) * stheta;
  h(2) = sin(phi) * stheta;
  h = h * H;
  //if(h(1)>0.001&&h(2)>=0){phi=atan(h(2)/h(1));}
  //if(h(1)>0.001&&h(2)<0){phi=atan(h(2)/h(1))+2*PI;}
  //if(h(1)<-0.001){phi=atan(h(2)/h(1))+PI;}
  jjj.mcalc(m,T,Hxc,h,jjj.Icalc_parstorage);
            //save physical properties of HT-point
    fprintf(fout,"%6.3f  %6.3f  %6.3f  %6.3f   %6.3f %6.3f %6.3f   %6.6f   %6.3f %6.3f %6.3f %6.3f\n",
           phi*180/PI,theta*180/PI,T,H,h(1),h(2),h(3),Norm(m),m(1),m(2),m(3),m*h/Norm(h));  
 }
}else{  // !do sipf
// as class par load  parameters from file
 if(verbose==1){printf("reading parameters from file mcphas.j\n");}
 char prefix [MAXNOFCHARINLINE];prefix[0]='\0';
 inipar ini("mcphas.ini",prefix); par inputpars("./mcphas.j"); inputpars.save("./results/_mcphas.j"); 
 nofthreads = ini.nofthreads;
  Vector Imax(1,inputpars.nofatoms*inputpars.nofcomponents);
  Vector Imom(1,inputpars.nofcomponents);
  Vector h1(1,inputpars.nofcomponents),h1ext(1,3);h1ext=0; 
  // here save single ion property files to results
  inputpars.save_sipfs("./results/_");
  //determine saturation momentum (used for scaling the plots, generation of qvectors)
  for(l=1;l<=inputpars.nofatoms;++l){h1=0;(*inputpars.jjj[l]).Icalc_parameter_storage_init(h1,h1ext,T); // initialize eigenstate matrix
  for (im=1;im<=inputpars.nofcomponents;++im){h1ext=0;h1=0;h1(im)=20*MU_B; //just put some high field
                            (*inputpars.jjj[l]).Icalc(Imom,T,h1,h1ext,z,u,(*inputpars.jjj[l]).Icalc_parstorage);
                            Imax(inputpars.nofcomponents*(l-1)+im)=Imom(im);
                                              }
                                  }
 // load testspinconfigurations (nooftstspinconfigurations,init-file,sav-file)
   testspincf testspins (ini.maxnoftestspincf,"./mcphas.tst","./results/mcphas.phs",inputpars.nofatoms,inputpars.nofcomponents);
   testspins.save("./results/_mcphas.tst","w");
   qvectors testqs (ini,inputpars.rez,Imax,"./results/mcphas.qvc",inputpars.nofatoms,inputpars.nofcomponents,verbose);
 // declare variable physprop (typa class physproperties)
   physproperties physprop(ini.nofspincorrs,ini.maxnofhkls,inputpars.nofatoms,inputpars.nofcomponents);
   	int nofstapoints=0,noffailedpoints=0,s=0;
 // loop different H /T points in phase diagram
/*
 for(double az=0;az<2*PI-0.00001;az+=2*PI/nofsteps)
 {h=H*(cos(az)*r1+sin(az)*r2);
  double phi=PI/2;if(h(2)<0){phi=3*PI/2;}
  if(h(1)>0.001&&h(2)>=0){phi=atan(h(2)/h(1));}
  if(h(1)>0.001&&h(2)<0){phi=atan(h(2)/h(1))+2*PI;}
  if(h(1)<-0.001){phi=atan(h(2)/h(1))+PI;}
  double theta=acos(h(3)/H);
*/
 for(int j=0; j<nofsteps; j++) {
  h(3) = 1. - (2. * j / nofsteps);
  double theta = acos(h(3));
  double phi = fmod(2. * j * PI_PHIINV, TWOPI);
  double stheta = sin(theta);
  h(1) = cos(phi) * stheta;
  h(2) = sin(phi) * stheta;
  h = h * H;
  //if(h(1)>0.001&&h(2)>=0){phi=atan(h(2)/h(1));}
  //if(h(1)>0.001&&h(2)<0){phi=atan(h(2)/h(1))+2*PI;}
  //if(h(1)<-0.001){phi=atan(h(2)/h(1))+PI;}
   // set field        
      physprop.T=T;
      physprop.H=h;
 //calculate physical properties at HT- point
   s=htcalc(physprop.H,T,ini,inputpars,testqs,testspins,physprop);
   if(s==1)break;
   //save physical properties of HT-point
   if(s==0)++nofstapoints;
   if(s==2)++noffailedpoints;
    fprintf(fout,"%6.3f  %6.3f  %6.3f  %6.3f   %6.3f %6.3f %6.3f   %6.6f   %6.3f %6.3f %6.3f %6.3f\n",
           phi*180/PI,theta*180/PI,T,H,h(1),h(2),h(3),Norm(physprop.m),physprop.m(1),physprop.m(2),physprop.m(3),physprop.m*h/Norm(h));  
  } // H/T loop 
   std::cout << "#!nofHTpoints=" << nofstapoints << "H-T points  successfully calculated" << std::endl;
   std::cout << "#!noffailedpoints=" << noffailedpoints << "H-T points in phasediagram failed to converge " << std::endl;
   std::cout << "#!fecalc - free energy calculation was attempted noffecalccalls=" << nofcalls << "times"  << std::endl;
   std::cout << "#!fecalc - free energy calculation was successful at noffecalcsuccess=" << successrate << "times"  << std::endl;
   std::cout << "#!fecalc - free energy diverged maxnofloopsDIV=" << nofmaxloopDIV << " times because maxnofloops was reached" << std::endl;
   std::cout << "#!fecalc - free energy diverged maxspinchangeDIV=" << nofmaxspinchangeDIV << " times because maxspinchange was reached" << std::endl;


} // do_sipffile
fclose(fout);
//  testspins.save(filemode);testqs.save(filemode);
   printf("RESULTS saved in directory ./results/  - files: anisotropy.out\n");
   double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
   std::cout << "#! Finished in cputime=" << cpu_duration << " seconds [CPU Clock] " << std::endl;
   
#ifdef _THREADS
std::cout << "#! nofthreads= " << nofthreads << " threads were used in parallel processing " << std::endl;
for (int ithread=0; ithread<nofthreads; ithread++) delete tin[ithread];
#else
std::cout << "# anisotropy was compiled without parallel processing option " << std::endl;
#endif

   fprintf(stderr,"**********************************************\n");
   fprintf(stderr,"          End of Program anisotropy\n");
   fprintf(stderr,"**********************************************\n");

return(0);
}

