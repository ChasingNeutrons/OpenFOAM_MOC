//Print areas for calculations by hand
OFstream myFile("integrals.txt");
	
forAll(entries, zoneI)
{
	const word& name = entries[zoneI].keyword();

	label zoneId = mesh.cellZones().findZoneID(name);
	if(zoneId == -1)
	{
	  Info << "Mesh zone: " << name << " does not exists" << endl;
	}

	myFile<<"Mesh zone: "<<name<<endl;
	myFile<<nl<<endl;
	scalar fluxInt=0;
	scalar fissInt=0;
	forAll(flux, energyI)
	{
		forAll(mesh.cellZones()[zoneId], cellIlocal)
		{

			label cellIglobal = mesh.cellZones()[zoneId][cellIlocal];
			
			fluxInt+=flux[energyI][cellIglobal]*area[cellIglobal];
			fissInt+=nuSigmaEff[energyI][cellIglobal]*flux[energyI][cellIglobal]*area[cellIglobal];
				
		}
		myFile<<"Energy "<<energyI<<nl<<endl;
		myFile<<fluxInt<<nl<<endl;
		fluxInt=0;
	}
	myFile<<"Fission source"<<nl<<endl;
	myFile<<fissInt<<nl<<endl;
			
}

