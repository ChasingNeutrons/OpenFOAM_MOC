//Fundamental mode rebalance to accelerate energy group convergence
//Occurs following an outer iteration to normalise flux and ensure neutron conservation.
//Performed on equivalent homogeneous system using flux and volume weighted XS
//Resulting group flux values then used to rebalance magnitude of heterogeneous fluxes

//ASSUMED ZERO LEAKAGE

//NEEDS CORRECTION TO HANDLE BOTH ANISOTROPIC AND LINEAR SOURCES

//First must generate homogeneous cross-sections
List<scalar> fluxFM(energyGroups);
List<scalar> prevFluxFM(energyGroups);
List<scalar> nuSigmaEffFM(energyGroups);
List<scalar> chiFM(energyGroups);
List<scalar> sigmaRFM(energyGroups);
List<scalar> qFM(energyGroups);
List<scalar> sigmaAFM(energyGroups);
List<List<scalar> > sigmaSFM(energyGroups);
forAll(sigmaSFM, energyI)
{
	sigmaSFM[energyI]=List<scalar>(energyGroups);
}

scalar fluxWeight;
scalar kfm;

//Normalisation factor to ensure only one neutron is absorbed in the system
scalar totalAbs;

//Calculate each cross section through flux weighting
forAll(flux, energyI)
{
	//Take initial guess for FM as the area weighted group fluxes
	fluxFM[energyI]=gSum(flux[energyI].internalField()*area);

	nuSigmaEffFM[energyI]=
	gSum(flux[energyI].internalField()*area*nuSigmaEff[energyI].internalField())/fluxFM[energyI];

	chiFM[energyI]=
	gSum(flux[energyI].internalField()*area*chi[energyI].internalField())/fluxFM[energyI];

	sigmaAFM[energyI]=
	gSum(flux[energyI].internalField()*area*sigmaA[energyI].internalField())/fluxFM[energyI];

	chiFM[energyI]=
	gSum(flux[energyI].internalField()*area*chi[energyI].internalField())/fluxFM[energyI];

	//sigmaR is the removal cross-section given by sigmaT - sigmaS g-g
	sigmaRFM[energyI]=
	gSum(flux[energyI].internalField()*area*
	(sigmaT[energyI].internalField()-sigmaS[energyI][energyI][0].internalField()))/fluxFM[energyI];


	//Should sigmaSg'-g be weighted by flux g' or g?
	forAll(flux, energyJ)
	{
		fluxWeight=gSum(flux[energyJ].internalField()*area);
		sigmaSFM[energyI][energyJ]=
		gSum(flux[energyJ].internalField()*sigmaS[energyI][energyJ][0].internalField()*area)/fluxWeight;
	}

	//Set previous fluxFM
	prevFluxFM[energyI]=fluxFM[energyI];
}

//Begin iterative process to calculate FM
scalar tolFM;
scalar testFM;
kfm=1.0;
do
{
	//Calculate and apply absorption normalisation
	totalAbs=0.0;
	forAll(flux, energyI)
	{
		totalAbs+=fluxFM[energyI]*sigmaAFM[energyI];
	}
	forAll(flux, energyI)
	{
		fluxFM[energyI]/=totalAbs;
	}

	//Calculate fundamental mode source
	forAll(flux,energyI)
	{
		//First calculate homogeneous source including fission contribution	
		qFM[energyI]=0.0;//chiFM[energyI];

		//Calculate contribution due to in-scatter (not including g-g scatter)
		forAll(flux,energyJ)
		{
			qFM[energyI]+=chiFM[energyI]*nuSigmaEffFM[energyJ]*fluxFM[energyJ]/kfm;
			if(energyI==energyJ){continue;}
			qFM[energyI]+=sigmaSFM[energyI][energyJ]*fluxFM[energyJ];
		}
	}

	//Update fluxFM
	forAll(fluxFM, energyI)
	{
		fluxFM[energyI]=qFM[energyI]/sigmaRFM[energyI];
	}

	//Check convergence
	tolFM=0.0;
	forAll(fluxFM, energyI)
	{
		testFM=mag(fluxFM[energyI]-prevFluxFM[energyI])/fluxFM[energyI];
		if(testFM>tolFM){tolFM=testFM;}
		prevFluxFM[energyI]=fluxFM[energyI];
	}

}while(tolFM>1e-6);

//Perform rebalance
forAll(fluxFM, energyI)
{
	Info<<fluxFM[energyI]/gSum(flux[energyI].internalField()*area)<<endl;
	flux[energyI]*=fluxFM[energyI]/gSum(flux[energyI].internalField()*area);

	if(anisotropy>0)
	{
		for(label l=0; l<anisotropy+1; l++)
		{
			for(label r=-l; r<l+1; r++)
			{
				fluxMo[energyI][l][r+l]*=fluxFM[energyI]/gSum(flux[energyI].internalField()*area);
			}
		}
	}
}






