forAll(fissionSource, energyI)
{
	fissionSource[energyI]*=0.0;
	forAll(fissionSource, energyJ)
	{
		fissionSource[energyI]+=chi[energyI] * nuSigmaEff[energyJ] * flux[energyJ];
	}
}
