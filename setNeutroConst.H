Info<<nl<<"Setting zone properties" << endl;
forAll(flux, energyI)
{
	forAll(entries, zoneI)
	{

		label zone = zoneI;

		const word& name = entries[zoneI].keyword();

		label zoneId = mesh.cellZones().findZoneID(name);
		if(zoneId == -1)
		{
		  Info << "Mesh zone: " << name << " does not exists" << endl;
		}


		forAll(mesh.cellZones()[zoneId], cellIlocal)
		{

			label cellIglobal = mesh.cellZones()[zoneId][cellIlocal];
			zone = zoneI;

			chi[energyI][cellIglobal] = chiList[zone][energyI];
			nuSigmaEff[energyI][cellIglobal] = nuSigmaEffList[zone][energyI];
			sigmaT[energyI][cellIglobal] = sigmaTList[zone][energyI];
			sigmaA[energyI][cellIglobal] = sigmaAList[zone][energyI];

			forAll(sigmaS[energyI],energyJ)
			{
				forAll(sigmaS[energyI][energyJ], l)
				{
					sigmaS[energyI][energyJ][l][cellIglobal] = sigmaFromToList[l][zone][energyJ][energyI];
				}
			}
		}
	}
}


