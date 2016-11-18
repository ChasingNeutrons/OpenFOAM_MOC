fissionSource[energyI] *= 0.0;
forAll(flux,energyJ)
{
	fissionSource[energyI] += chi[energyI] * nuSigmaEff[energyJ] * flux[energyJ];
}
