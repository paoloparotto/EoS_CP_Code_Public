/* Copyright (c) 2018, Paolo Parotto and Debora Mroczek, 
 * Department of Physics, University of Houston, Houston, TX 77204, US. */

/* This main produces an EoS matching Lattice QCD at muB=0, and containing a critical 
 * point in the 3D Ising universality class, in a parametrized form.
 * It allows for different choices of constraints on the shape of the critical 
 * line, which reduce the number of parameters. */ 

#define NRANSI

// --------- INCLUDE -------- // 
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_deriv.h>
#include "../include/nrD.h"
#include "../include/nrutilD.h"
#include "../include/Variables.h"
#include "../include/Functions.h"

// --------- DEFINE --------- //
#define PI 3.141592653589
#define N 2


// The main
int main(int argc, char *argv[])
{
	char buff[FILENAME_MAX];

	/* Set phase diagram region */
	lowT=5; highT=821;
	lowMU=0; highMU=601;
	lowT_out=30; highT_out=800;
	lowMU_out=0; highMU_out=450;

// --------- ALLOCATE ------------------- // 
{{{
	// Vectors for the root finding. 
	x=vector(1,N);	f=vector(1,N);
	// Matrix for Jacobian. 
	JJ=matrix(1,2,1,2);
	
	// Rank 3 tensor for coordinates.
	Coords=f3tensor(lowT,highT,0,2*highMU,1,2);
	// Vectors for Chi's at muB=0. 
	Chi0LatVec=vector(lowT,highT);			
	Chi2LatVec=vector(lowT,highT);			
	Chi4LatVec=vector(lowT,highT);			
	Chi0IsingVec=vector(lowT,highT);			
	Chi2IsingVec=vector(lowT,highT);			
	Chi4IsingVec=vector(lowT,highT);			
	Chi0NoIsingVec=vector(lowT,highT);		
	Chi2NoIsingVec=vector(lowT,highT);		
	Chi4NoIsingVec=vector(lowT,highT);		
	dChi0NoIsingdTVec=vector(lowT,highT);	
	dChi2NoIsingdTVec=vector(lowT,highT);	
	dChi4NoIsingdTVec=vector(lowT,highT);	
	d2Chi0NoIsingdT2Vec=vector(lowT,highT);	
	d2Chi2NoIsingdT2Vec=vector(lowT,highT);	
	d2Chi4NoIsingdT2Vec=vector(lowT,highT);	
	
	// Matrices for thermodynamic functions over the phase diagram. 
	// NonIsing
	PressNoIsingMat=matrix(lowT,highT,0,highMU);			
	dPressNoIsingdTMat=matrix(lowT,highT,0,highMU);				
	dPressNoIsingdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressNoIsingdT2Mat=matrix(lowT,highT,0,highMU);		
	d2PressNoIsingdmuB2Mat=matrix(lowT,highT,0,highMU);			
	d2PressNoIsingdTdmuBMat=matrix(lowT,highT,0,highMU);
	// NonIsing (after filtering)
	PressNoIsingFilterMat=matrix(lowT,highT,0,highMU);		
	dPressNoIsingFilterdTMat=matrix(lowT,highT,0,highMU);		
	dPressNoIsingFilterdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressNoIsingFilterdT2Mat=matrix(lowT,highT,0,highMU);	
	d2PressNoIsingFilterdmuB2Mat=matrix(lowT,highT,0,highMU);	
	d2PressNoIsingFilterdTdmuBMat=matrix(lowT,highT,0,highMU);
	d3PressNoIsingFilterdmuB3Mat=matrix(lowT,highT,0,highMU); 
	d4PressNoIsingFilterdmuB4Mat=matrix(lowT,highT,0,highMU);
	// Ising
	PressIsingMat=matrix(lowT,highT,0,highMU);				
	dPressIsingdTMat=matrix(lowT,highT,0,highMU);				
	dPressIsingdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressIsingdT2Mat=matrix(lowT,highT,0,highMU);			
	d2PressIsingdmuB2Mat=matrix(lowT,highT,0,highMU);			
	d2PressIsingdTdmuBMat=matrix(lowT,highT,0,highMU);
	d3PressIsingdmuB3Mat=matrix(lowT,highT,0,highMU);		
	d4PressIsingdmuB4Mat=matrix(lowT,highT,0,highMU);
	// Ising + NonIsing (after filtering)
	PressTotMat=matrix(lowT,highT,0,highMU);				
	dPressTotdTMat=matrix(lowT,highT,0,highMU);					
	dPressTotdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressTotdT2Mat=matrix(lowT,highT,0,highMU);			
	d2PressTotdmuB2Mat=matrix(lowT,highT,0,highMU);				
	d2PressTotdTdmuBMat=matrix(lowT,highT,0,highMU);
	d3PressTotdmuB3Mat=matrix(lowT,highT,0,highMU);			
	d4PressTotdmuB4Mat=matrix(lowT,highT,0,highMU);	
	// HRG
	PressHRGMat=matrix(lowT,highT,0,highMU);				
	dPressHRGdTMat=matrix(lowT,highT,0,highMU);					
	dPressHRGdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressHRGdT2Mat=matrix(lowT,highT,0,highMU);			
	d2PressHRGdmuB2Mat=matrix(lowT,highT,0,highMU);				
	d2PressHRGdTdmuBMat=matrix(lowT,highT,0,highMU);
	d3PressHRGdmuB3Mat=matrix(lowT,highT,0,highMU);			
	d4PressHRGdmuB4Mat=matrix(lowT,highT,0,highMU);	
	// Total + HRG
	PressTotHRGMat=matrix(lowT,highT,0,highMU);				
	dPressTotHRGdTMat=matrix(lowT,highT,0,highMU);				
	dPressTotHRGdmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressTotHRGdT2Mat=matrix(lowT,highT,0,highMU);		
	d2PressTotHRGdmuB2Mat=matrix(lowT,highT,0,highMU);			
	d2PressTotHRGdTdmuBMat=matrix(lowT,highT,0,highMU);
	d3PressTotHRGdmuB3Mat=matrix(lowT,highT,0,highMU);		
	d4PressTotHRGdmuB4Mat=matrix(lowT,highT,0,highMU);
	// Final (normalized)
	PressFinalMat=matrix(lowT,highT,0,highMU);				
	EntropyFinalMat=matrix(lowT,highT,0,highMU);				
	BarDensityFinalMat=matrix(lowT,highT,0,highMU);
	EnerDensityFinalMat=matrix(lowT,highT,0,highMU);		
	SpSoundFinalMat=matrix(lowT,highT,0,highMU);				
	Chi2FinalMat=matrix(lowT,highT,0,highMU);			
	Chi4FinalMat=matrix(lowT,highT,0,highMU);
	// If LAT only is chosen (not normalized)
	PressLATonlyMat=matrix(lowT,highT,0,highMU);			
	dPressLATonlydTMat=matrix(lowT,highT,0,highMU);				
	dPressLATonlydmuBMat=matrix(lowT,highT,0,highMU);	
	d2PressLATonlydT2Mat=matrix(lowT,highT,0,highMU);		
	d2PressLATonlydmuB2Mat=matrix(lowT,highT,0,highMU);			
	d2PressLATonlydTdmuBMat=matrix(lowT,highT,0,highMU);
	d3PressLATonlydmuB3Mat=matrix(lowT,highT,0,highMU);		
	d4PressLATonlydmuB4Mat=matrix(lowT,highT,0,highMU);			
	PressLATonlyFilterMat=matrix(lowT,highT,0,highMU);
	// If LAT only is chosen (normalized)
	PressLATonlyNormMat=matrix(lowT,highT,0,highMU);		
	EntropyLATonlyNormMat=matrix(lowT,highT,0,highMU);			
	BarDensityLATonlyNormMat=matrix(lowT,highT,0,highMU);	
	EnerDensityLATonlyNormMat=matrix(lowT,highT,0,highMU);	
	SpSoundLATonlyNormMat=matrix(lowT,highT,0,highMU);			
	Chi2LATonlyNormMat=matrix(lowT,highT,0,highMU);
	Chi4LATonlyNormMat=matrix(lowT,highT,0,highMU);	
}}}

 	/* Assign the name of the main folder where the program lives and 
	 * the files we wish to import are located. */
	getcwd(buff,FILENAME_MAX);
	printf("Current working directory is: \n\n%s \n\n",buff);

// --------- IMPORT LATTICE CHIS -------- //
{{{	
	// Chi's from Lattice are imported and stored in a vector.
	printf("Importing lattice data...\n");
	FILE *LatIn = fopen("input/Lattice_Data_5_821_dT1.dat","r");
	if (LatIn == 0){
  	fprintf(stderr, "failed to open Lattice Data\n");
  	exit(1);
 	}
	else printf("File imported succesfully!\n");

 	// Save lattice chis into vectors (make them dimension of energy^4). Print to file latout if desired.	
 	//	FILE *Latout = fopen("Chis_Lat_5_821_dT1.dat","w");
	for(i=lowT;fscanf(LatIn,"%lf %lf %lf %lf",&xIn1,&xIn2,&xIn3,&xIn4) !=EOF;i++){
  	Chi0LatVec[i] = xIn2*pow(i,4);
   	Chi2LatVec[i] = xIn3*pow(i,4);
   	Chi4LatVec[i] = xIn4*pow(i,4);
   	//fprintf(Latout,"%3.1f %12.16f	%12.16f	%12.16f\n",(double) i,Chi0LatVec[i],Chi2LatVec[i],Chi4LatVec[i]);
 	}
 	fclose(LatIn);	//fclose(Latout);
 	printf("Successfully stored lattice data.\n");
}}}

// --------- IMPORT PARAMETERS -----------//
{{{	
	// Open the parameter file
 	FILE *ParametersIn = fopen(argv[1], "r");
 	if (ParametersIn == 0){
 		fprintf(stderr,"failed to open paremeters file\n");
 		exit(1);
 	}
	else printf("File opened correctly.\n");

	// Import and store parameters
 	for(i=5;fscanf(ParametersIn,"%s %lf %lf %lf %lf %lf %lf", v0, &v1, &v2, &v3, &v4, &v5, &v6) !=EOF; i++){
 		const char * MODE = v0;
 		const char * PATH = p0;

 		printf("MODE: %s\n", MODE);

 		if(strcmp(MODE,"LAT") == 0){
 			strcpy(MODESTR, "LAT");
 			TC = v1;    muBC = v2;  angle1 = v3;    angle2 = v4;    ww = v5;    rho = v6;
 			printf("You chose to insert no critical point; the Taylor expansion from Lattice QCD will be returned.\n");
 		} else if(strcmp(MODE,"FREE") == 0){
 			strcpy(MODESTR, "FREE");
 			TC = v1;    muBC = v2;  angle1 = v3;    angle2 = v4;    ww = v5;    rho = v6;
 			printf("The parameters you entered are:\nTC = %f\nmuBC = %f\nangle1 = %f\nangle2 = %f\nw = %f\nrho = %f\n",TC,muBC,angle1,angle2,ww,rho);
 		} else if(strcmp(MODE,"PAR") == 0){
 			strcpy(MODESTR, "PAR");
 			T0 = v1;    kappa = v2; muBC = v3;  anglediff = v4; ww = v5;    rho = v6;
 			TC = T0 + kappa/T0 * muBC * muBC; angle1 = 180/PI*fabs(atan(-2.0*kappa/T0*muBC)); angle2 = angle1 + anglediff;
 			printf("PARABOLA\n");
 			printf("The parameters you entered are:\nTC = %f\nmuBC = %f\nangle1 = %f\nangle2 = %f\nw = %f\nrho = %f\n",TC,muBC,angle1,angle2,ww,rho);
 		} else if(strcmp(MODE,"STR") == 0){
 			strcpy(MODESTR, "STR");
 			T0 = v1;    muBC = v2;  angle1 = v3;    angle2 = v4;    ww = v5;    rho = v6;
 			TC = T0 - atan(angle1*PI/180)*muBC;
 			printf("The parameters you entered are:\nTC = %f\nmuBC = %f\nangle1 = %f\nangle2 = %f\nw = %f\nrho = %f\n",TC,muBC,angle1,angle2,ww,rho);
 		} else {};
     	break;
 	}
	// Close parameter file.
 	fclose(ParametersIn);
}}}
	printf("Finished importing and setting parameters.\n\n");

// ----------------------- START LATTICE CASE --------------------//
{{{ 
	/* Output Taylor expansion from Lattice QCD without critical point if MODE = LAT.*/
  if(strcmp(MODESTR,"LAT") == 0){

	printf("You chose to consider the lattice-only case, without a Critical Point.\n");
  // Creating directory and moving into it.
	mkdir("LATonly", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	chdir("LATonly");

	// Calculate Taylor expanded pressure.
	for (i=lowT;i<=highT; i++) for(j=0; j<=highMU; j++){
		muBval = (double) j; 	Tval = (double) i;
  	PressLATonlyMat[i][j] = Chi0LatVec[i] + 1.0/2.0*Chi2LatVec[i]*pow(muBval/Tval,2) + 1.0/24.0*Chi4LatVec[i]*pow(muBval/Tval,4);
	}

	// Take all the needed derivatives of the pressure.
	// Wrt T
	deriv_matrix(PressLATonlyMat,dPressLATonlydTMat,1,lowT,highT,0,highMU,0);	
	deriv_matrix(dPressLATonlydTMat,d2PressLATonlydT2Mat,1,lowT,highT,0,highMU,0);
	// Wrt muB
	deriv_matrix(PressLATonlyMat,dPressLATonlydmuBMat,2,lowT,highT,0,highMU,0);
	deriv_matrix(dPressLATonlydmuBMat,d2PressLATonlydmuB2Mat,2,lowT,highT,0,highMU,1);
	deriv_matrix(d2PressLATonlydmuB2Mat,d3PressLATonlydmuB3Mat,2,lowT,highT,0,highMU,0);
	deriv_matrix(d3PressLATonlydmuB3Mat,d4PressLATonlydmuB4Mat,2,lowT,highT,0,highMU,1);
	// Wrt T and muB
	deriv_matrix(dPressLATonlydTMat,d2PressLATonlydTdmuBMat,2,lowT,highT,0,highMU,0);
	
	// That is it. Now we just need to combine these into the final observables (normalized).
	printf("Calculating thermodynamics quantities. \n");

	// Create the files for export.
	FILE *FilePressLATonlyNorm = fopen("PressLATonlyNorm3D.dat", "w");		FILE *FileEntrLATonlyNorm = fopen("EntrLATonlyNorm3D.dat", "w");
	FILE *FileBarDensLATonlyNorm = fopen("BarDensLATonlyNorm3D.dat", "w");	FILE *FileEnerDensLATonlyNorm = fopen("EnerDensLATonlyNorm3D.dat", "w");
	FILE *FileSpSoundLATonlyNorm = fopen("SpsoundLATonlyNorm3D.dat", "w");	FILE *FileChi2LATonlyNorm = fopen("Chi2LATonlyNorm3D.dat", "w");
	FILE *FileChi4LATonlyNorm = fopen("Chi4LATonlyNorm3D.dat", "w");
		
	// Export.
	for (i=lowT_out; i<=highT_out; i++) for(j=0;j<=highMU_out; j++){
		Tval = (double) i; muBval = (double) j;
		
		PressLATonlyNormMat[i][j] = PressLATonlyMat[i][j]/pow(Tval,4);
 		EntropyLATonlyNormMat[i][j] = dPressLATonlydTMat[i][j]/pow(Tval,3);
  	BarDensityLATonlyNormMat[i][j] = (dPressLATonlydmuBMat[i][j]/pow(Tval,3))*cure;
  	EnerDensityLATonlyNormMat[i][j] = (dPressLATonlydTMat[i][j]*Tval - PressLATonlyMat[i][j] + muBval*dPressLATonlydmuBMat[i][j])/pow(Tval,4);
  	SpSoundLATonlyNormMat[i][j] = (dPressLATonlydmuBMat[i][j]*dPressLATonlydmuBMat[i][j]*d2PressLATonlydT2Mat[i][j] - 2.0*dPressLATonlydTMat[i][j]
  											*dPressLATonlydmuBMat[i][j]*d2PressLATonlydTdmuBMat[i][j] + dPressLATonlydTMat[i][j]
   											*dPressLATonlydTMat[i][j]*d2PressLATonlydmuB2Mat[i][j])
   											*1.0/(dPressLATonlydTMat[i][j]*Tval + muBval*dPressLATonlydmuBMat[i][j])
   											*1.0/(d2PressLATonlydT2Mat[i][j]*d2PressLATonlydmuB2Mat[i][j]-d2PressLATonlydTdmuBMat[i][j]
  											*d2PressLATonlydTdmuBMat[i][j]);     								
   	Chi2LATonlyNormMat[i][j] = d2PressLATonlydmuB2Mat[i][j]/pow(Tval,2);
   	Chi4LATonlyNormMat[i][j] = d4PressLATonlydmuB4Mat[i][j];

  	fprintf(FilePressLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, PressLATonlyNormMat[i][j]);
  	fprintf(FileEntrLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, EntropyLATonlyNormMat[i][j]);
  	fprintf(FileBarDensLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, BarDensityLATonlyNormMat[i][j]);
  	fprintf(FileEnerDensLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, EnerDensityLATonlyNormMat[i][j]);
 	  fprintf(FileSpSoundLATonlyNorm,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, SpSoundLATonlyNormMat[i][j]);
		fprintf(FileChi2LATonlyNorm,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, Chi2LATonlyNormMat[i][j]);
		fprintf(FileChi4LATonlyNorm,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, Chi4LATonlyNormMat[i][j]);
  	}
	// Close the files.
	fclose(FilePressLATonlyNorm);		fclose(FileEntrLATonlyNorm);		fclose(FileBarDensLATonlyNorm);
	fclose(FileEnerDensLATonlyNorm);	fclose(FileSpSoundLATonlyNorm);		fclose(FileChi2LATonlyNorm);	
	fclose(FileChi4LATonlyNorm);

	printf("The procedure is completed. Exiting.\n");

	return 0;
  }
}}}
// ----------------------- END LATTICE CASE --------------------- //


// ----------------------- START GENERAL CASE ------------------- //

	// We can immediately determine dTC and dmuBC from TC, muBC, ww and rho. 
	dTC = ww*TC; 
	dmuBC = dTC*rho;

	/* Determine the filenames according to the parameter choice made. */
	{{{
	sprintf(nameFolder, "Files_%s_%d_%d_%d_%d_%d_%d",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameCoords, "Coords_%s_%d_%d_%d_%d_%d_%d.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChisIsing, "Chis_Ising_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChisNoIsing, "Chis_No_Ising_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedChisNoIsingdT, "dChis_No_Ising_dT_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2ChisNoIsingdT2, "d2Chis_No_Ising_dT2%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressIsing3D, "Press_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTIsing3D, "dPdT_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBIsing3D, "dPdmuB_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2Ising3D, "d2PdT2_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2Ising3D, "d2PdmuB2_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named4PdmuB4Ising3D, "d4PdmuB4_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBIsing3D, "d2PdTdmuB_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressNoIsing3D, "Press_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTNoIsing3D, "dPdT_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBNoIsing3D, "dPdmuB_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2NoIsing3D, "d2PdT2_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2NoIsing3D, "d2PdmuB2_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBNoIsing3D, "d2PdTdmuB_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressIsingPlusNoIsing3D, "Press_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTIsingPlusNoIsing3D, "dPdT_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBIsingPlusNoIsing3D, "dPdmuB_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2IsingPlusNoIsing3D, "d2PdT2_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2IsingPlusNoIsing3D, "d2PdmuB2_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBIsingPlusNoIsing3D, "d2PdTdmuB_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressTotHRG3D, "Press_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTTotHRG3D, "dPdT_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBTotHRG3D, "dPdmuB_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2TotHRG3D, "d2PdT2_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2TotHRG3D, "d2PdmuB2_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBTotHRG3D, "d2PdTdmuB_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressFinal3D, "Press_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameEntrFinal3D, "Entr_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameBarDensFinal3D, "BarDens_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameEnerDensFinal3D, "EnerDens_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameSpsoundFinal3D, "SpSound_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChi2Final3D, "Chi2_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChi4Final3D, "Chi4_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	}}}

	// Move into folder for output
	chdir("output");
	// Create folder with name depending on parameters, and move to the folder.
	mkdir(nameFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	chdir(nameFolder);

	// For calculations hereafter, angles are expressed in radians. 
	angle1 = angle1*PI/180;
	angle2 = angle2*PI/180;

	// Once all the parameters are set, the Jacobian matrix of the (r,h) -> (T,muB) map is generated and
	// a value is obtained for drdmuB and dhdmuB.  
	Jacobian(JJ,dTC,dmuBC,angle1,angle2);
	drdmuB = JJ[1][2];
	dhdmuB = JJ[2][2];
	drdT = JJ[1][1];
	dhdT = JJ[2][1];
	printf("\nThe Jacobian calculation gave:\ndrdmuB = %12.16f\ndhdT = %12.16f\ndrdT = %12.16f\ndhmuB = %12.16f",drdmuB,dhdmuB,drdT,dhdT);

// ----------- GENERATE COORDINATES -------- //
{{{
	// Now everything is set, the program generates the coordinate file, assigning it the name stored in 
	//	the string "namecoords".  
	FILE *FileCoords=fopen(nameCoords,"w");
	printf("\n\nCreating the coordinates file...\n");
	// Time is taken at start and end to measure how long it took to generate the coordinates file.
	time(&start);
	for(Ti=lowT;Ti<=highT;Ti+=1.0) for(muBi=-highMU;muBi<=highMU;muBi+=1.0){
		h_i = hh(Ti,muBi);
		r_i  = rr(Ti, muBi);

    /* Condition for choosing init condition in root finding. */
   	x[1] = 1.0;
		if(h_i > 0) x[2] = 1.0;
		if(h_i <= 0) x[2] = -1.0;
						
		/* The actual root finding routine. */
		newt(x,N,&check,funcv);
		funcv(N,x,f);
		if (check) printf("Convergence problems.\n");
		/* Print to file the coordinates {T,muB,R,Theta} (in this order). */
		fprintf(FileCoords,"%3.1f %3.1f %2.16f %2.16f\n",Ti,muBi,x[1],x[2]);
	} 
	time(&end);
	/* A control on whether the loops reached the end, and did not get stuck, otherwise file is eliminated. */	
	if(Ti == (highT+1) && muBi == (highMU+1))	printf("\nThe file %s was created successfully in %d seconds\n", nameCoords, (int) difftime(end,start));
	else{
		printf("\nThe file %s was NOT created successfully, removing it..\n",nameCoords);
		remove(nameCoords);
	} 								
	fclose(FileCoords);

	// Coordinates are imported and stored in a Rank 3 tensor Coords. 
	printf("\nImporting coordinates from file %s \n",nameCoords);
	FileCoords=fopen(nameCoords,"r");
	for(j=0, x1int=lowT-1;fscanf(FileCoords,"%lf %lf %lf %lf\n", &xIn1, &xIn2, &xIn3, &xIn4) != EOF; j++ ){
		x2int = j % (2*highMU+1);
		if (x2int == 0) x1int++;
 		Coords[x1int][x2int][1] = xIn3;
		Coords[x1int][x2int][2] = xIn4; 
	}
	fclose(FileCoords);
  printf("\nImported coordinates successfully.\n");
}}}

// ----------- ISING ----------------------- //
{{{
	/* The coordinates corresponding to muB=0 are used to generate the critical Chi's at different temperatures, and exported. */
	printf("\nGenerating Ising at muB=0\n");
	FILE *FileChiIsing=fopen(nameChisIsing,"w");
	for(i=lowT;i<=highT;i++){
		Tval = (double) i;
		Chi0IsingVec[i] = -G(Coords[i][highMU][1],Coords[i][highMU][2])*pow(TC,4);
		Chi2IsingVec[i] = -Tval*Tval*d2GdmuB2ConT(Coords[i][highMU][1],Coords[i][highMU][2])*pow(TC,4);
		Chi4IsingVec[i] = -Tval*Tval*Tval*Tval*d4GdmuB4ConT(Coords[i][highMU][1],Coords[i][highMU][2])*pow(TC,4);

		fprintf(FileChiIsing,"%3.1f %12.16f %12.16f %12.16f\n",(double) i,Chi0IsingVec[i],Chi2IsingVec[i],Chi4IsingVec[i]);
	}
	fclose(FileChiIsing);
	printf("\nGeneration successful\n");

	/* The pressure is generated over the whole range of coordinates.
	 * NOTE: the pressure is symmetrized around muB=0 to ensure all odd order 
	 * derivatives vanish at muB=0. 
	 * NOTE: file output is currently disabled. (Un)comment file opening and 
	 * closing, as well as the export line in the loop, to turn on/off the 
	 * file export.*/
	printf("\nGenerating Ising in 3D\n");
	//FILE *FilePressIsing=fopen(namePressIsing3D,"w");				FILE *FiledPdTIsing=fopen(namedPdTIsing3D,"w");
	//FILE *FiledPdmuBIsing=fopen(namedPdmuBIsing3D,"w");			FILE *Filed2PdT2Ising=fopen(named2PdT2Ising3D,"w");
	//FILE *Filed2PdmuB2Ising=fopen(named2PdmuB2Ising3D,"w");		FILE *Filed2PdTdmuBIsing=fopen(named2PdTdmuBIsing3D,"w");
	//FILE *Filed4PdmuB4Ising=fopen(named4PdmuB4Ising3D,"w");
	for(i=lowT;i<highT;i++) for(j=0;j<=highMU;j++){
		k = j + highMU;
		if(j==0){
			PressIsingMat[i][j] = - G(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			dPressIsingdTMat[i][j] = - dGdTConmuB(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			dPressIsingdmuBMat[i][j] = 0.0;
			d2PressIsingdT2Mat[i][j] = - d2GdT2ConmuB(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			d2PressIsingdmuB2Mat[i][j] = - d2GdmuB2ConT(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			d2PressIsingdTdmuBMat[i][j] = 0.0;
			//d3PressIsingdmuB3Mat[i][j] = - d3GdmuB3ConT(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			//d4PressIsingdmuB4Mat[i][j] = - d4GdmuB4ConT(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
		} else{
	 		PressIsingMat[i][j] = - (G(Coords[i][k][1],Coords[i][k][2]) + G(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		dPressIsingdTMat[i][j] = - (dGdTConmuB(Coords[i][k][1],Coords[i][k][2]) + dGdTConmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		dPressIsingdmuBMat[i][j] = - (dGdmuBConT(Coords[i][k][1],Coords[i][k][2]) - dGdmuBConT(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		d2PressIsingdT2Mat[i][j] = - (d2GdT2ConmuB(Coords[i][k][1],Coords[i][k][2]) + d2GdT2ConmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		d2PressIsingdmuB2Mat[i][j] = - (d2GdmuB2ConT(Coords[i][k][1],Coords[i][k][2]) + d2GdmuB2ConT(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		d2PressIsingdTdmuBMat[i][j] = - (d2GdTdmuB(Coords[i][k][1],Coords[i][k][2]) - d2GdTdmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		//d3PressIsingdmuB3Mat[i][j] = - (d3GdmuB3ConT(Coords[i][k][1],Coords[i][k][2]) - d3GdmuB3ConT(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);	
			//d4PressIsingdmuB4Mat[i][j] = - (d4GdmuB4ConT(Coords[i][k][1],Coords[i][k][2]) + d4GdmuB4ConT(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);	 			 			
 		}
			
 		//fprintf(FilePressIsing,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, PressIsingMat[i][j]);
 		//fprintf(FiledPdTIsing,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, dPressIsingdTMat[i][j]);
 		//fprintf(FiledPdmuBIsing,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, dPressIsingdmuBMat[i][j]);
 		//fprintf(Filed2PdT2Ising,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, d2PressIsingdT2Mat[i][j]);
		//fprintf(Filed2PdmuB2Ising,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, d2PressIsingdmuB2Mat[i][j]);
 		//fprintf(Filed2PdTdmuBIsing,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, d2PressIsingdTdmuBMat[i][j]);
 		//fprintf(Filed4PdmuB4Ising,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, d4PressIsingdmuB4Mat[i][j]);	 		
	}
	//fclose(FilePressIsing); fclose(FiledPdTIsing); fclose(FiledPdmuBIsing); 
	//fclose(Filed2PdT2Ising); fclose(Filed2PdmuB2Ising);	fclose(Filed2PdTdmuBIsing);		
	//fclose(Filed4PdmuB4Ising);
}}}

	// Remove coordinate file
	remove(nameCoords);

// ----------- NON-ISING	------------------ //
{{{
	/* Non-Ising Chi's are calculated, exported and stored. 
	 * NOTE: file output is currently disabled. (Un)comment file opening and 
	 * closing, as well as the export line in the loop,	to turn on/off the 
	 * file export. */
	printf("\nGenerating non-Ising at muB = 0\n");
	//FILE *FileChiNoIsing = fopen(nameChisNoIsing,"w");
	for (i=lowT; i<=highT; i++){
	    Tval = (double) i;
	    Chi0NoIsingVec[i] = Chi0LatVec[i] - Chi0IsingVec[i];
	    Chi2NoIsingVec[i] = Chi2LatVec[i] - Chi2IsingVec[i];
	    Chi4NoIsingVec[i] = Chi4LatVec[i] - Chi4IsingVec[i];
        
	   	//fprintf(FileChiNoIsing,"%3.1f %12.16f %12.16f %12.16f\n", Tval, Chi0NoIsingVec[i],Chi2NoIsingVec[i],Chi4NoIsingVec[i]);
	}
	//fclose(FileChiNoIsing);

	/* The non-Ising pressure in 3D is calculated from Taylor expansion, 
	 * stored and exported.  
	 * NOTE: file output is currently disabled. (Un)comment file opening and 
	 * closing, as well as the export line in the loop, to turn on/off the 
	 * file export.*/
 	printf("\nCalculating non-Ising in 3D.\n");
	//FILE *FilePressNoIsing=fopen(namePressNoIsing3D,"w");			
	for (i=lowT; i<=highT; i++) for(j=0; j<=highMU; j++){
		muBval = (double) j; Tval = (double) i;
    PressNoIsingMat[i][j] = Chi0NoIsingVec[i]+(1.0/2.0)*Chi2NoIsingVec[i]*pow(muBval/Tval,2)
												+(1.0/24.0)*Chi4NoIsingVec[i]*pow(muBval/Tval,4);

		//fprintf(FilePressNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,PressNoIsingMat[i][j]);
	}
	//fclose(FilePressNoIsing); 

	/* Filtering of regular pressure here. Widths of the filter are also set here. 
	 * Result is stored and exported. 
	 * NOTE: file output is currently disabled. (Un)comment file opening and 
	 * closing, as well as the export line in the loop,	to turn on/off the 
	 * file export.*/
	printf("\nFiltering non-Ising Pressure in 3D\n");
	sigmax = 18.0;	sigmay = 1.0;
	//FILE *FilePressNoIsingFilter = fopen("NoIsingPressFilter.dat","w");
	for (i=lowT; i<=highT; i++)	for(j=0; j<=highMU; j++){
		muBval = (double) j; Tval = (double) i;

		kmin = i-sigmax; kmax = i+sigmax;
		lmin = j-sigmay; lmax = j+sigmay;	
 
		sum = 0.0; norm = 0.0;
		for(k=kmin;k<=kmax;k++) for(l=lmin;l<=lmax;l++){
			if(k<=lowT)		kint = lowT;
			else if(k>=highT)	kint = highT;	
			else		 	kint = k;
					
			if(l<=0)  		lint = -l;
			else if(l>=highMU)	lint = highMU;		
			else			lint = l;
									
			sum+=GausFunc(kint-i,lint-j,sigmax,sigmay)*PressNoIsingMat[kint][lint];
			norm+=GausFunc(kint-i,lint-j,sigmax,sigmay);
		}
		PressNoIsingFilterMat[i][j] = sum/norm;
		//fprintf(FilePressNoIsingFilter, "%3.1f %3.1f %12.16f\n", muBval,Tval,PressNoIsingFilterMat[i][j]);
	}
	//fclose(FilePressNoIsingFilter);

	/* Derivatives (numerical) of the filtered Non-Ising pressure and all thermodynamic 
	 * quantities of interest are calculated, and stored, over the phase diagram. */
	printf("\nCalculating (numerical) derivatives of the filtered Non-Ising Pressure\n");
	//FILE *FiledPdTNoIsingFilter = fopen("NoIsingdPdTFilter.dat","w");						
	//FILE *FiledPdmuBNoIsingFilter = fopen("NoIsingdPdmuBFilter.dat","w");
	//FILE *Filed2PdT2NoIsingFilter = fopen("NoIsingd2PdT2Filter.dat","w");				
	//FILE *Filed2PdmuB2NoIsingFilter = fopen("NoIsingd2PdmuB2Filter.dat","w");
	//FILE *Filed2PdTdmuBNoIsingFilter = fopen("NoIsingd2PdTdmuBFilter.dat","w");	
	//FILE *Filed3PdmuB3NoIsingFilter = fopen("NoIsingd3PdmuB3Filter.dat","w");
	//FILE *Filed4PdmuB4NoIsingFilter = fopen("NoIsingd4PdmuB4Filter.dat","w");

	// Take all the needed derivatives
	// Wrt T
	deriv_matrix(PressNoIsingFilterMat,dPressNoIsingFilterdTMat,1,lowT,highT,0,highMU,0);	
	deriv_matrix(dPressNoIsingFilterdTMat,d2PressNoIsingFilterdT2Mat,1,lowT,highT,0,highMU,0);
	// Wrt muB
	deriv_matrix(PressNoIsingFilterMat,dPressNoIsingFilterdmuBMat,2,lowT,highT,0,highMU,0);
	deriv_matrix(dPressNoIsingFilterdmuBMat,d2PressNoIsingFilterdmuB2Mat,2,lowT,highT,0,highMU,1);
	deriv_matrix(d2PressNoIsingFilterdmuB2Mat,d3PressNoIsingFilterdmuB3Mat,2,lowT,highT,0,highMU,0);
	deriv_matrix(d3PressNoIsingFilterdmuB3Mat,d4PressNoIsingFilterdmuB4Mat,2,lowT,highT,0,highMU,1);
	// Wrt T and muB
	deriv_matrix(dPressNoIsingFilterdTMat,d2PressNoIsingFilterdTdmuBMat,2,lowT,highT,0,highMU,0);
}}}

// ----------- ISING + NON-ISING ----------- //
{{{
	/* Calculate pressure and derivatives summing Ising and non-Ising 
	 * pressure, then normalize. Store and export. 
	 * NOTE: file output is currently disabled. (Un)comment file opening and 
	 * closing, as well as the export line in the loop,	to turn on/off the 
	 * file export.*/
  printf("\nCalculating other derivatives in 3D \n");
  //FILE *FilePressTot = fopen(namePressIsingPlusNoIsing3D, "w");					
	//FILE *FiledPdTIsingPlusNoIsing=fopen(namedPdTIsingPlusNoIsing3D,"w");
	//FILE *FiledPdmuBIsingPlusNoIsing=fopen(namedPdmuBIsingPlusNoIsing3D,"w");		
	//FILE *Filed2PdT2IsingPlusNoIsing=fopen(named2PdT2IsingPlusNoIsing3D,"w");
	//FILE *Filed2PdmuB2IsingPlusNoIsing=fopen(named2PdmuB2IsingPlusNoIsing3D,"w");	
	//FILE *Filed2PdTdmuBIsingPlusNoIsing=fopen(named2PdTdmuBIsingPlusNoIsing3D,"w");
	for (i=lowT; i<=highT; i++) for(j=0;j<=highMU; j++){
		Tval = (double) i; muBval = (double) j;
		
		PressTotMat[i][j] = PressIsingMat[i][j] + PressNoIsingFilterMat[i][j];
   	dPressTotdTMat[i][j] = dPressIsingdTMat[i][j] + dPressNoIsingFilterdTMat[i][j];
   	dPressTotdmuBMat[i][j] = dPressIsingdmuBMat[i][j] + dPressNoIsingFilterdmuBMat[i][j];
   	d2PressTotdT2Mat[i][j] = d2PressIsingdT2Mat[i][j] + d2PressNoIsingFilterdT2Mat[i][j];
   	d2PressTotdmuB2Mat[i][j] = d2PressIsingdmuB2Mat[i][j] + d2PressNoIsingFilterdmuB2Mat[i][j];  		
   	d2PressTotdTdmuBMat[i][j] = d2PressIsingdTdmuBMat[i][j] + d2PressNoIsingFilterdTdmuBMat[i][j];   
   	d3PressTotdmuB3Mat[i][j] = d3PressIsingdmuB3Mat[i][j] + d3PressNoIsingFilterdmuB3Mat[i][j];  		
   	d4PressTotdmuB4Mat[i][j] = d4PressIsingdmuB4Mat[i][j] + d4PressNoIsingFilterdmuB4Mat[i][j];  		

   	//fprintf(FilePressTot,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, PressTotMat[i][j]);
 		//fprintf(FiledPdTIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,dPressTotdTMat[i][j]);
 		//fprintf(FiledPdmuBIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,dPressTotdmuBMat[i][j]);
 		//fprintf(Filed2PdT2IsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressTotdT2Mat[i][j]);
 		//fprintf(Filed2PdmuB2IsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressTotdmuB2Mat[i][j]);
 		//fprintf(Filed2PdTdmuBIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressTotdTdmuBMat[i][j]);
  }
	//fclose(FilePressTot);					
	//fclose(FiledPdTIsingPlusNoIsing); 		
	//fclose(FiledPdmuBIsingPlusNoIsing); 
	//fclose(Filed2PdT2IsingPlusNoIsing);	
	//fclose(Filed2PdmuB2IsingPlusNoIsing);	
	//fclose(Filed2PdTdmuBIsingPlusNoIsing);
}}}

// ----------- MERGING WITH HRG -------------//
{{{
 	/* Now, the smooth joining with the HRG pressure starts. Import HRG pressure and store. */
	printf("\nImporting HRG Pressure \n");
	FILE *FilePressHRG = fopen("../../input/Press_HRG_MUB000601_T005300_dT1.dat", "r");
	if (FilePressHRG == 0){
		fprintf(stderr, "failed to open HRG Pressure \n");
		exit(1);
	}
	for(i=0, x2int=0;fscanf(FilePressHRG,"%lf %lf %lf\n",&xIn1,&xIn2,&xIn3) !=EOF;i++){
		x1int = (i % 817) + 5;
		PressHRGMat[x1int][x2int] = xIn3*pow(x1int,4);
		if(x1int == 821) x2int++;
	}
	fclose(FilePressHRG); 

 	/* Calculate derivatives of the HRG pressure. Store and export. */
	/* NOTE: file output is currently disabled. (Un)comment file opening and 
	 * closing, as well as the export line in the loop,	to turn on/off the 
	 * file export.*/
  printf("\nCalculating derivatives of HRG pressure in 3D \n");
	// Take all the needed derivatives
	{{{
	//FILE *FiledPdTHRG=fopen("dPdT_HRG_3D.dat","w");
	//FILE *FiledPdmuBHRG=fopen("dPdmuB_HRG_3D.dat","w");		FILE *Filed2PdT2HRG=fopen("d2PdT2_HRG_3D.dat","w");
	//FILE *Filed2PdmuB2HRG=fopen("d2PdmuB2_HRG_3D.dat","w");	FILE *Filed2PdTdmuBHRG=fopen("d2PdTdmuB_HRG_3D.dat","w");
	//FILE *Filed3PdmuB3HRG=fopen("d3PdmuB3_HRG_3D.dat","w");	FILE *Filed4PdmuB4HRG=fopen("d4PdmuB4_HRG_3D.dat","w");
	// Wrt T
	deriv_matrix(PressHRGMat,dPressHRGdTMat,1,lowT,highT,0,highMU,0);	
	deriv_matrix(dPressHRGdTMat,d2PressHRGdT2Mat,1,lowT,highT,0,highMU,0);
	// Wrt muB
	deriv_matrix(PressHRGMat,dPressHRGdmuBMat,2,lowT,highT,0,highMU,0);
	deriv_matrix(dPressHRGdmuBMat,d2PressHRGdmuB2Mat,2,lowT,highT,0,highMU,1);
	deriv_matrix(d2PressHRGdmuB2Mat,d3PressHRGdmuB3Mat,2,lowT,highT,0,highMU,0);
	deriv_matrix(d3PressHRGdmuB3Mat,d4PressHRGdmuB4Mat,2,lowT,highT,0,highMU,1);
	// Wrt T and muB
	deriv_matrix(dPressHRGdTMat,d2PressHRGdTdmuBMat,2,lowT,highT,0,highMU,0);
	}}}

	/* Final pressure, joined with the HRG one is calculated, stored and exported. 
	 * The parameter deltaT for merging with HRG is also set here. */
	/* NOTE: file output is currently disabled. (Un)comment file opening and 
	 * closing, as well as the export line in the loop,	to turn on/off the 
	 * file export.*/
	printf("\nMerging with HRG data.\n");
  /* Merging with HRG data. */
	{{{
	//FILE *FilePressTotHRG = fopen(namePressTotHRG3D, "w");		FILE *FiledPdTTotHRG=fopen(namedPdTTotHRG3D,"w");
	//FILE *FiledPdmuBTotHRG=fopen(namedPdmuBTotHRG3D,"w");		FILE *Filed2PdT2TotHRG=fopen(named2PdT2TotHRG3D,"w");
	//FILE *Filed2PdmuB2TotHRG=fopen(named2PdmuB2TotHRG3D,"w");	FILE *Filed2PdTdmuBTotHRG=fopen(named2PdTdmuBTotHRG3D,"w");

	// Set delta T
	deltaT = 17.0;

	/* The actual merging */
	for (i=lowT; i<=highT; i++) for(j=0;j<=highMU; j++){
		Tval = (double) i; muBval = (double) j;
  	Targum = (Tval - (T0 + kappa/T0*muBval*muBval - 23.0))/deltaT;

  	if(Targum >= 4.5)    		PressTotHRGMat[i][j] = PressTotMat[i][j];
  	else if(Targum <= -4.5) PressTotHRGMat[i][j] = PressHRGMat[i][j];
  	else                    PressTotHRGMat[i][j] = PressTotMat[i][j]*0.5*(1.0 + tanh(Targum)) 
                                  + PressHRGMat[i][j]*0.5*(1.0 - tanh(Targum));
	
  	if(Targum >= 4.5)       dPressTotHRGdTMat[i][j] = dPressTotdTMat[i][j];
  	if(Targum <= -4.5)      dPressTotHRGdTMat[i][j] = dPressHRGdTMat[i][j];
  	else                    dPressTotHRGdTMat[i][j] = dPressTotdTMat[i][j]*0.5*(1.0 + tanh(Targum))
                            			+ dPressHRGdTMat[i][j]*0.5*(1.0 - tanh(Targum))
                                 	+ PressTotMat[i][j]*0.5/deltaT*pow(cosh(Targum),-2)
                                 	- PressHRGMat[i][j]*0.5/deltaT*pow(cosh(Targum),-2);       

  	if(Targum >= 4.5)       dPressTotHRGdmuBMat[i][j] = dPressTotdmuBMat[i][j];
  	if(Targum <= -4.5)      dPressTotHRGdmuBMat[i][j] = dPressHRGdmuBMat[i][j];
  	else                    dPressTotHRGdmuBMat[i][j] = dPressTotdmuBMat[i][j]*0.5*(1.0 + tanh(Targum))
                            			+ dPressHRGdmuBMat[i][j]*0.5*(1.0 - tanh(Targum))
                                 	- PressTotMat[i][j]*0.5*pow(cosh(Targum),-2)*(2*kappa/T0*muBval)/deltaT
                                 	+ PressHRGMat[i][j]*0.5*pow(cosh(Targum),-2)*(2*kappa/T0*muBval)/deltaT;

  	if(Targum >= 4.5)       d2PressTotHRGdT2Mat[i][j] = d2PressTotdT2Mat[i][j];
  	if(Targum <= -4.5)      d2PressTotHRGdT2Mat[i][j] = d2PressHRGdT2Mat[i][j];
  	else                    d2PressTotHRGdT2Mat[i][j] = d2PressTotdT2Mat[i][j]*0.5*(1.0 + tanh(Targum)) 
                            			+ d2PressHRGdT2Mat[i][j]*0.5*(1.0 - tanh(Targum))
                           				+ dPressTotdTMat[i][j]/deltaT*pow(cosh(Targum),-2)
                                 	- dPressHRGdTMat[i][j]/deltaT*pow(cosh(Targum),-2)
                                 	- PressTotMat[i][j]*pow(cosh(Targum),-2)*tanh(Targum)/pow(deltaT,2) 
                                 	+ PressHRGMat[i][j]*pow(cosh(Targum),-2)*tanh(Targum)/pow(deltaT,2);   

  	if(Targum >= 4.5)       d2PressTotHRGdmuB2Mat[i][j] = d2PressTotdmuB2Mat[i][j];
  	if(Targum <= -4.5)      d2PressTotHRGdmuB2Mat[i][j] = d2PressHRGdmuB2Mat[i][j];
   	else                    d2PressTotHRGdmuB2Mat[i][j] = d2PressTotdmuB2Mat[i][j]*0.5*(1.0 + tanh(Targum))
                            			+ d2PressHRGdmuB2Mat[i][j]*0.5*(1.0 - tanh(Targum))
                            			- dPressTotdmuBMat[i][j]*pow(cosh(Targum),-2)*(2*kappa/T0*muBval)/deltaT
                           				+ dPressHRGdmuBMat[i][j]*pow(cosh(Targum),-2)*(2*kappa/T0*muBval)/deltaT
                                 	- PressTotMat[i][j]*0.5*pow(cosh(Targum),-2)*(2*tanh(Targum)
                                 		*(2*kappa/T0*muBval)/deltaT
                                 		*(2*kappa/T0*muBval)/deltaT
                                 		+ 2*kappa/T0/deltaT)
                                 	+ PressHRGMat[i][j]*0.5*pow(cosh(Targum),-2)*(2*tanh(Targum)
                                 		*(2*kappa/T0*muBval)/deltaT
                                 		*(2*kappa/T0*muBval)/deltaT
                                 		+ 2*kappa/T0/deltaT);

  	if(Targum >= 4.5)       d2PressTotHRGdTdmuBMat[i][j] = d2PressTotdTdmuBMat[i][j];
  	if(Targum <= -4.5)      d2PressTotHRGdTdmuBMat[i][j] = d2PressHRGdTdmuBMat[i][j];
  	else                    d2PressTotHRGdTdmuBMat[i][j] = d2PressTotdTdmuBMat[i][j]*0.5*(1.0 + tanh(Targum))
                            			+ d2PressHRGdTdmuBMat[i][j]*0.5*(1.0 - tanh(Targum))
                                 	+ dPressTotdmuBMat[i][j]*0.5/deltaT*pow(cosh(Targum),-2)
                                 	- dPressHRGdmuBMat[i][j]*0.5/deltaT*pow(cosh(Targum),-2) 
                            			- dPressTotdTMat[i][j]*0.5/deltaT*(2*kappa/T0*muBval)*pow(cosh(Targum),-2)
                                 	+ dPressHRGdTMat[i][j]*0.5/deltaT*(2*kappa/T0*muBval)*pow(cosh(Targum),-2)
                                 	+ PressTotMat[i][j]*pow(cosh(Targum),-2)/pow(deltaT,2)
                                 		*tanh(Targum)*(2*kappa/T0*muBval)
                                 	- PressHRGMat[i][j]*pow(cosh(Targum),-2)/pow(deltaT,2)
                                 		*tanh(Targum)*(2*kappa/T0*muBval);    

		// Tanh merging is not done when considering the 4th order derivative
		d4PressTotHRGdmuB4Mat[i][j] = d4PressTotdmuB4Mat[i][j];
    
		/* Chi4 HRG merging */
		{{{	
		/*		if(Targum >= 4.5)	 	d4PressTotHRGdmuB4Mat[i][j] = d4PressTotdmuB4Mat[i][j];
        	if(Targum <= -4.5)		d4PressTotHRGdmuB4Mat[i][j] = d4PressHRGdmuB4Mat[i][j];
  				else 					d4PressTotHRGdmuB4Mat[i][j] = d4PressTotdmuB4Mat[i][j]*0.5*(1.0 + tanh(Targum))
            								+ d4PressHRGdmuB4Mat[i][j]*0.5*(1.0 - tanh(Targum))
            								- d3PressTotdmuB3Mat[i][j]*4.0/deltaT*(kappa/T0*muBval)*pow(cosh(Targum),-2) 
            								+ d3PressHRGdmuB3Mat[i][j]*4.0/deltaT*(kappa/T0*muBval)*pow(cosh(Targum),-2)
            								- d2PressTotdmuB2Mat[i][j]*6.0*kappa/(deltaT*T0)*pow(cosh(Targum),-2)
            									*(1.0 + 4.0*kappa/(deltaT*T0)*muBval*muBval*tanh(Targum))
            								+ d2PressHRGdmuB2Mat[i][j]*6.0*kappa/(deltaT*T0)*pow(cosh(Targum),-2)
            									*(1.0 + 4.0*kappa/(deltaT*T0)*muBval*muBval*tanh(Targum))
            								- dPressTotdmuBMat[i][j]*16.0*kappa/(deltaT*T0)*kappa/(deltaT*T0)
            									*muBval*pow(cosh(Targum),-2)*(- 2.0*kappa/(deltaT*T0)*muBval*muBval*pow(cosh(Targum),-2)
            										+ 3.0*tanh(Targum) + 4.0*kappa/(deltaT*T0)*muBval*muBval*pow(tanh(Targum),2))  
            								+ dPressHRGdmuBMat[i][j]*16.0*kappa/(deltaT*T0)*kappa/(deltaT*T0)
            									*muBval*pow(cosh(Targum),-2)*(- 2.0*kappa/(deltaT*T0)*muBval*muBval*pow(cosh(Targum),-2)
            										+ 3.0*tanh(Targum) + 4.0*kappa/(deltaT*T0)*muBval*muBval*pow(tanh(Targum),2))
            								- PressTotMat[i][j]*kappa/(deltaT*T0)*kappa/(deltaT*T0)*pow(cosh(Targum),-2)
            									*(- 48.0*kappa/(deltaT*T0)*pow(muBval,2)*pow(cosh(Targum),-2) 
            										+ 12.0*tanh(Targum) - 128.0*pow(kappa/(deltaT*T0),2)*pow(muBval,4)*pow(cosh(Targum),-2)*tanh(Targum)
            										+ 96.0*kappa/(deltaT*T0)*pow(muBval,2)*pow(tanh(Targum),2)
            										+ 64.0*pow(kappa/(deltaT*T0),2)*pow(muBval,4)*pow(tanh(Targum),3))   
            								+ PressHRGMat[i][j]*kappa/(deltaT*T0)*kappa/(deltaT*T0)*pow(cosh(Targum),-2)
            									*(- 48.0*kappa/(deltaT*T0)*pow(muBval,2)*pow(cosh(Targum),-2) 
            										+ 12.0*tanh(Targum) - 128.0*pow(kappa/(deltaT*T0),2)*pow(muBval,4)*pow(cosh(Targum),-2)*tanh(Targum)
            										+ 96.0*kappa/(deltaT*T0)*pow(muBval,2)*pow(tanh(Targum),2)
            										+ 64.0*pow(kappa/(deltaT*T0),2)*pow(muBval,4)*pow(tanh(Targum),3));       */                 			
		}}}

		/* Output is disabled now */
		{{{
   	//fprintf(FilePressTotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, PressTotHRGMat[i][j]);
  	//fprintf(FiledPdTTotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, dPressTotHRGdTMat[i][j]);
   	//fprintf(FiledPdmuBTotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, dPressTotHRGdmuBMat[i][j]);
   	//fprintf(Filed2PdT2TotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, d2PressTotHRGdT2Mat[i][j]);
 		//fprintf(Filed2PdmuB2TotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, d2PressTotHRGdmuB2Mat[i][j]);
  	//fprintf(Filed2PdTdmuBTotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, d2PressTotHRGdTdmuBMat[i][j]);
		}}}            
	}
	//fclose(FilePressTotHRG);	fclose(FiledPdTTotHRG); 		fclose(FiledPdmuBTotHRG); 
	//fclose(Filed2PdT2TotHRG);	fclose(Filed2PdmuB2TotHRG);		fclose(Filed2PdTdmuBTotHRG);
  }}}
}}}
	printf("\nMerged with HRG.\n");

// ----------- OUTPUT ---------------------- //
{{{
	/* Now all the derivatives have been merged with the HRG, we can finish with 
	 * the thermodynamics and export the quantities we want, in particular we 
	 * will normalize with suitable powers of the temperature to obtain 
	 * dimensionless quantities.*/
	printf("\nMerging with HRG complete. \n \nCalculating thermodynamics quantities. \n");
	/* Calculate thermodynamic quantities */
	{{{
	FILE *FilePressFinal = fopen(namePressFinal3D, "w");			FILE *FileEntrFinal = fopen(nameEntrFinal3D, "w");
	FILE *FileBarDensFinal = fopen(nameBarDensFinal3D, "w");	FILE *FileEnerDensFinal = fopen(nameEnerDensFinal3D, "w");
	FILE *FileSpSoundFinal = fopen(nameSpsoundFinal3D, "w");	FILE *FileChi2Final = fopen(nameChi2Final3D, "w");
	FILE *FileChi4Final = fopen(nameChi4Final3D, "w");
	for (i=lowT_out; i<=highT_out; i++) for(j=0;j<=highMU_out; j++){
		Tval = (double) i; muBval = (double) j;

   	PressFinalMat[i][j] = PressTotHRGMat[i][j]/pow(Tval,4);
   	EntropyFinalMat[i][j] = dPressTotHRGdTMat[i][j]/pow(Tval,3);
   	BarDensityFinalMat[i][j] = dPressTotHRGdmuBMat[i][j]/pow(Tval,3);
   	EnerDensityFinalMat[i][j] = (dPressTotHRGdTMat[i][j]*Tval - PressTotHRGMat[i][j] 
																+ muBval*dPressTotHRGdmuBMat[i][j])/pow(Tval,4);
   	SpSoundFinalMat[i][j] = (dPressTotHRGdmuBMat[i][j]*dPressTotHRGdmuBMat[i][j]*d2PressTotHRGdT2Mat[i][j] 
																- 2.0*dPressTotHRGdTMat[i][j]*dPressTotHRGdmuBMat[i][j]*d2PressTotHRGdTdmuBMat[i][j] 
																+ dPressTotHRGdTMat[i][j]*dPressTotHRGdTMat[i][j]*d2PressTotHRGdmuB2Mat[i][j])
   													*1.0/(dPressTotHRGdTMat[i][j]*Tval + muBval*dPressTotHRGdmuBMat[i][j])
   													*1.0/(d2PressTotHRGdT2Mat[i][j]*d2PressTotHRGdmuB2Mat[i][j]
																-d2PressTotHRGdTdmuBMat[i][j]*d2PressTotHRGdTdmuBMat[i][j]);
    Chi2FinalMat[i][j] = d2PressTotHRGdmuB2Mat[i][j]/pow(Tval,2);
    Chi4FinalMat[i][j] = d4PressTotHRGdmuB4Mat[i][j];

		fprintf(FilePressFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, PressFinalMat[i][j]);
 		fprintf(FileEntrFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, EntropyFinalMat[i][j]);
 		fprintf(FileBarDensFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, BarDensityFinalMat[i][j]);
 		fprintf(FileEnerDensFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, EnerDensityFinalMat[i][j]);
	 	fprintf(FileSpSoundFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, SpSoundFinalMat[i][j]);
		fprintf(FileChi2Final,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, Chi2FinalMat[i][j]);
 		fprintf(FileChi4Final,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, Chi4FinalMat[i][j]);
 	}
	fclose(FilePressFinal);			fclose(FileEntrFinal);			fclose(FileBarDensFinal);
	fclose(FileEnerDensFinal);	fclose(FileSpSoundFinal);		
	fclose(FileChi2Final);			fclose(FileChi4Final);
	}}}
}}}	
	printf("\nProcedure completed.\n");
// ------------------------ END GENERAL CASE ---------------------- //


// ----------- FREE ALLOCATIONS ------------ //
	printf("\nFree memory allocation.\n\n");
{{{
	/* Remove coordinate file */
	remove(nameCoords);
  /* Vectors for root finding */
	free_vector(x,1,N);
	free_vector(f,1,N);
  /* Jacobian */
	free_matrix(JJ,1,2,1,2);
  /* Tensor for coordinates*/
	free_f3tensor(Coords,lowT,highT,0,2*highMU,1,2);
  /* Vectors */
	free_vector(Chi0LatVec,lowT,highT);			
	free_vector(Chi2LatVec,lowT,highT);			
	free_vector(Chi4LatVec,lowT,highT);		
	free_vector(Chi0IsingVec,lowT,highT);		
	free_vector(Chi2IsingVec,lowT,highT);		
	free_vector(Chi4IsingVec,lowT,highT);	
	free_vector(Chi0NoIsingVec,lowT,highT);	
	free_vector(Chi2NoIsingVec,lowT,highT);	
	free_vector(Chi4NoIsingVec,lowT,highT);	
	/* For LAT-only mode */ 
	free_matrix(PressLATonlyMat,lowT,highT,0,highMU);					free_matrix(PressLATonlyFilterMat,lowT,highT,0,highMU);		free_matrix(dPressLATonlydTMat,lowT,highT,0,highMU);		
	free_matrix(dPressLATonlydmuBMat,lowT,highT,0,highMU);		free_matrix(d2PressLATonlydT2Mat,lowT,highT,0,highMU);		free_matrix(d2PressLATonlydmuB2Mat,lowT,highT,0,highMU);	
	free_matrix(d2PressLATonlydTdmuBMat,lowT,highT,0,highMU);	free_matrix(d3PressLATonlydmuB3Mat,lowT,highT,0,highMU);	free_matrix(d4PressLATonlydmuB4Mat,lowT,highT,0,highMU);
  /* For LAT-only mode, normalized */
	free_matrix(PressLATonlyNormMat,lowT,highT,0,highMU);				free_matrix(EntropyLATonlyNormMat,lowT,highT,0,highMU);		free_matrix(BarDensityLATonlyNormMat,lowT,highT,0,highMU);
	free_matrix(EnerDensityLATonlyNormMat,lowT,highT,0,highMU);	free_matrix(SpSoundLATonlyNormMat,lowT,highT,0,highMU);		free_matrix(Chi2LATonlyNormMat,lowT,highT,0,highMU);
	free_matrix(Chi4LATonlyNormMat,lowT,highT,0,highMU);
  /* Ising */ 
	free_matrix(PressNoIsingMat,lowT,highT,0,highMU);				free_matrix(dPressNoIsingdTMat,lowT,highT,0,highMU);			free_matrix(dPressNoIsingdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressNoIsingdT2Mat,lowT,highT,0,highMU);	free_matrix(d2PressNoIsingdmuB2Mat,lowT,highT,0,highMU);	free_matrix(d2PressNoIsingdTdmuBMat,lowT,highT,0,highMU);
  /* No-Ising filtered */
	free_matrix(PressNoIsingFilterMat,lowT,highT,0,highMU);				free_matrix(dPressNoIsingFilterdTMat,lowT,highT,0,highMU);		
	free_matrix(dPressNoIsingFilterdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressNoIsingFilterdT2Mat,lowT,highT,0,highMU);	free_matrix(d2PressNoIsingFilterdmuB2Mat,lowT,highT,0,highMU);	free_matrix(d2PressNoIsingFilterdTdmuBMat,lowT,highT,0,highMU);
	free_matrix(d3PressNoIsingFilterdmuB3Mat,lowT,highT,0,highMU);	free_matrix(d4PressNoIsingFilterdmuB4Mat,lowT,highT,0,highMU);
  /* Isig */
	free_matrix(PressIsingMat,lowT,highT,0,highMU);			    free_matrix(dPressIsingdTMat,lowT,highT,0,highMU);			free_matrix(dPressIsingdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressIsingdT2Mat,lowT,highT,0,highMU);		free_matrix(d2PressIsingdmuB2Mat,lowT,highT,0,highMU);		free_matrix(d2PressIsingdTdmuBMat,lowT,highT,0,highMU);
	free_matrix(d3PressIsingdmuB3Mat,lowT,highT,0,highMU);	    free_matrix(d4PressIsingdmuB4Mat,lowT,highT,0,highMU);
  /* Total */
	free_matrix(PressTotMat,lowT,highT,0,highMU);			    free_matrix(dPressTotdTMat,lowT,highT,0,highMU);			free_matrix(dPressTotdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressTotdT2Mat,lowT,highT,0,highMU);		    free_matrix(d2PressTotdmuB2Mat,lowT,highT,0,highMU);		free_matrix(d2PressTotdTdmuBMat,lowT,highT,0,highMU);
	free_matrix(d3PressTotdmuB3Mat,lowT,highT,0,highMU);		free_matrix(d4PressTotdmuB4Mat,lowT,highT,0,highMU);
  /* HRG */
	free_matrix(PressHRGMat,lowT,highT,0,highMU);			    free_matrix(dPressHRGdTMat,lowT,highT,0,highMU);			free_matrix(dPressHRGdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressHRGdT2Mat,lowT,highT,0,highMU);		    free_matrix(d2PressHRGdmuB2Mat,lowT,highT,0,highMU);		free_matrix(d2PressHRGdTdmuBMat,lowT,highT,0,highMU);
	free_matrix(d3PressHRGdmuB3Mat,lowT,highT,0,highMU);		free_matrix(d4PressHRGdmuB4Mat,lowT,highT,0,highMU);
  /* Total+HRG */
	free_matrix(PressTotHRGMat,lowT,highT,0,highMU);		    free_matrix(dPressTotHRGdTMat,lowT,highT,0,highMU);			free_matrix(dPressTotHRGdmuBMat,lowT,highT,0,highMU);
	free_matrix(d2PressTotHRGdT2Mat,lowT,highT,0,highMU);	    free_matrix(d2PressTotHRGdmuB2Mat,lowT,highT,0,highMU);		free_matrix(d2PressTotHRGdTdmuBMat,lowT,highT,0,highMU);
	free_matrix(d3PressTotHRGdmuB3Mat,lowT,highT,0,highMU);	    free_matrix(d4PressTotHRGdmuB4Mat,lowT,highT,0,highMU);
	/* FInal output */
	free_matrix(PressFinalMat,lowT,highT,0,highMU);			    free_matrix(EntropyFinalMat,lowT,highT,0,highMU);			free_matrix(BarDensityFinalMat,lowT,highT,0,highMU);
	free_matrix(EnerDensityFinalMat,lowT,highT,0,highMU);	    free_matrix(SpSoundFinalMat,lowT,highT,0,highMU);			free_matrix(Chi2FinalMat,lowT,highT,0,highMU);	
	free_matrix(Chi4FinalMat,lowT,highT,0,highMU);
}}}

  /* End */
	return 0;
}


#undef NRANSI
