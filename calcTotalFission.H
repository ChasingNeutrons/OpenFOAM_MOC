totalFissions*=0.0;
forAll(flux,energyJ)
{
	totalFissions += nuSigmaEff[energyJ] * flux[energyJ];
}

