/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    MOClinearAni

Description

    2D MOC solver for rectangular geometry. Solves multi-group, multi-region problems.
    (Soon) Capable of dealing with anisotropic and linear source problems. 

    Given a cuboidal mesh, traces rays from the bottom y-boundary and left and right
    x-boundaries. Ray tracing information is stored and used during the transport sweep
    to calculate group scalar fluxes and criticality.

    Settings for solver are specified in MOCsettings file in the system folder:
    npo = number of polar angles (must be between 1 and 3)
    naz = number of azimuthal angles (input will be adjusted to be divisible by 4)
    spacing = unadjusted spacing between tracks in metres (OpenFOAM standard unit)
    tolerance = maximum allowable residual before solver terminates

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//Included for indexedOctree ray tracing
#include "polyMesh.H"
#include "meshTools.H"
#include "treeDataFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Define constants

const scalar PI=constant::mathematical::pi;
const scalar one_over_4_PI=1.0/(4.0*PI);
#include "sphericalHarmonic.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readSettings.H"
    #include "readNuclearDataExt.H"
    #include "createFields.H"
    #include "setNeutroConst.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    //Create polar and azimuthal information
    Info<< "Creating quadrature\n" << endl;
    #include "initialiseQuad.H"

    //Initialising ray information storage
    Info<< "Creating ray tracing storage\n" << endl;
    #include "createRayStorage.H"
    
    //Ray Tracing using indexedOctree
    Info<<"Tracing rays across geometry\n" << endl;
    #include "rayTrace.H" 

    //ADD THIS
    //Output ray start and end information if so desired
    //Info<<"Outputting ray start and end points"<<endl;
    //#include "outputRayCoords.H"
    
    //Matching rays for applying BCs
    Info<<"Matching complementary rays\n" << endl;
    #include "rayMatching.H"

    //Define other parameters used
    scalar delta, albedo, afluxIn, wt, wa, sina, cosa, sinT, sigT, A, tau, err2, e1, e2, e3;
    scalar oldK=1.0;
    scalar err=1.0;
    scalar qb, qh, mu;
    label segNum, ind, i0;
    int c1, c2, m, refRay, refDir;

    //Doesn't work! Need to investigate!
    //Homogeneous acceleration
    //Info<<"Initialising flux using energy rebalance"<<endl; 
    //#include "FMRebalance.H"

    //Initialise fission source, flux, angular flux and
    //previous iterate of each
    #include "initialiseFluxSource.H"

    //ADD THIS!!!
    //Correct volumes and areas to account for unstructured mesh
    //if a radius is set in readSettings
    //Info<<"Correcting fuel volume & area for unstructured mesh"<<endl;
    //#include "geometryCorrect.H" 

    //Compute and store exponential factors 
    Info<<"Pre-computing exponential quantities"<<endl;
    #include "computeExp.H"

    //Alternative to compute exponentials - ADD THIS
    //Info<<"Tabulating exponential lookup table"<<endl;
    //#include "tabulateExponentialTable.H"

    if(linearSource)
    {
	Info<<"Pre-computing C and M factors"<<endl;
	if(anisotropy==0)
	{		    
		//Compute the C factors for linear source calculations
		#include "computeCFactors.H"
    	}else{
		//Compute anisotropic C factors for lienar source
		#include "computeAniCFactors.H"
	}
    }
    
    //Begin loop
    Info<< "\nStarting iterations\n" << endl;
    while (runTime.loop() && err>tol)
    {
	Info<< "Iteration " << runTime.timeName() << endl;
	err=0.0;        

	//Transport sweep - branches for either linear or flat source and isotropic or anisotropic scattering
	if(linearSource)
	{
		if(anisotropy==0)
		{
			//Update source moments
			#include "calcSourceMoments.H"	
	
			//TRANSPORT SWEEP
			#include "transportSweepLinearSource.H"	
		}else{
			//Update anisotropic source moments
			#include "calcAnisotropicSourceMoments.H"	
	
			//TRANSPORT SWEEP with anisotropic source
			#include "transportSweepLinearAngularFlux.H"

			//Update scalar flux moments
			#include "fluxCalculateLinear.H"
		}		
	}else{
		if(anisotropy==0)
		{
			//TRANSPORT SWEEP
			#include "transportSweep.H"
		}else{
			//Transport Sweep with anisotropic source
			#include "transportSweepAngularFlux.H"

			//Update scalar flux moments
			#include "fluxCalculate.H"
		}			
	}
	
	//Doesn't work! Need to investigate!
	//Homogeneous energy group acceleration
	//#include "FMRebalance.H"

	//ADD THIS!!!
	//Chebyshev acceleration
	//#include "chebyshevAcceleration.H"

	//Calculate new total fission source
	#include "calcTotalFission.H"
	volFissions=gSum(totalFissions.internalField() * area);

        //Update k-eff using power iterations
	oldK=keff;
	keff *= volFissions/prevFissions;
	prevFissions=volFissions;
	
        //Q update and calculate error
	if(anisotropy==0)
	{
		forAll(Q,energyI)
		{
			#include "calcFissionSource.H"	
			#include "calcScatteringSource.H"
			Q[energyI]= scatteringSource[energyI] + fissionSource[energyI]/keff;
		}
	}else{
		#include "calcAnisotropicSource.H"
		#include "calcFissionSourceAni.H"
	}

	//Calculate error - on fission source rather than flux perhaps?
	e2=0.0;
	scalar eT=0.0;
	forAll(flux, energyI)
	{
		eT=max(mag(flux[energyI].internalField()-oldFlux[energyI].internalField())/flux[energyI].internalField());
		if(eT>e2){e2=eT;}
		oldFlux[energyI]=flux[energyI];
	}
	e3=mag(keff-oldK)/mag(keff);
	if(e2>e3){err=e2;}
	else {err=e3;}

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"<<nl
            << "K-eff = " << keff << nl
            << "Error in k = " << e3 << nl
	    << "Max error in scalar flux = "<< e2 << nl << endl;
    }  

    if (!runTime.loop())
    {
	Info<<"Calculation did not converge"<<endl;
    }
   
    //For multi-group problems, print ratio of group fluxes
    if(energyGroups>1)
    {
	scalar fluxRatio;
	forAll(flux, energyI)
	{
		for(int energyJ=energyI+1; energyJ<energyGroups; energyJ++)
		{
			fluxRatio=gSum(flux[energyJ].internalField()*vol)/gSum(flux[energyI].internalField()*vol);
			Info<<"Group"<<energyJ+1<<"/Group"<<energyI+1<<" = "<<fluxRatio<<endl;
		}
	}
    }

    //Include a total fission source rather than energy specific!
    forAll(flux, energyI)
    {
	flux[energyI].write();
	if(energyI>0){
		fissionSource[0]+=fissionSource[energyI];
	}
    }
    fissionSource[0].write();

    //Integrate group fluxes in each mesh zone
    if (zoneIntegral)
    {
	#include "integrateZoneFlux.H"
    }
	
   
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
