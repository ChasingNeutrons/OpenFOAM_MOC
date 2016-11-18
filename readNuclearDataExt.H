	Info<< "\nReading nuclear data\n" << endl;

	IOdictionary nuclearData
	(
		IOobject
		(
			"nuclearData",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	label energyGroups(nuclearData.lookupOrDefault("energyGroups",1));
	label anisotropy(nuclearData.lookupOrDefault("anisotropy",0)); //cant presently go higher than 4 due to Legendre function
	scalar keff(nuclearData.lookupOrDefault("keff",1.0));

	PtrList<entry> entries(nuclearData.lookup("zones"));
	label zoneNumber = entries.size();

	//Initialise nuclear data 
	PtrList<scalarField > chiList(zoneNumber);
	PtrList<scalarField > nuSigmaEffList(zoneNumber);
	PtrList<scalarField > sigmaTList(zoneNumber);
	PtrList<scalarField > sigmaAList(zoneNumber);
	//Define scattering matrices for order of anisotropy
	PtrList<PtrList<scalarSquareMatrix> > sigmaFromToList(anisotropy+1);
	for(label l=0; l<anisotropy+1; l++)			
	{
		sigmaFromToList.set(l, new PtrList<scalarSquareMatrix>(zoneNumber));	
	}

        //Fill the data
        forAll(entries,zoneI)
        {
                dictionary& dict = entries[zoneI].dict();
		
		const word& name = entries[zoneI].keyword();
		word scatterOrder;

		Info << "name: " << name << endl;

		//Nuclear data
		chiList.set(zoneI,new scalarField("chi",dict,energyGroups));
		nuSigmaEffList.set(zoneI,new scalarField("nuSigmaEff",dict,energyGroups));
		sigmaTList.set(zoneI,new scalarField("sigmaT",dict,energyGroups));
		sigmaAList.set(zoneI,new scalarField("sigmaA",dict,energyGroups));
		forAll(sigmaFromToList, l)
		{
			//Read the scattering matrix of appropriate order
			if(l==0){scatterOrder="";}
			else{scatterOrder=Foam::name(l);}

			sigmaFromToList[l].set(zoneI,new scalarSquareMatrix(dict.lookup("scatteringMatrix"+scatterOrder)));
		}
        }

	Info << "Finished reading nuclear data " << endl;
