scatteringSource[energyI]*=0.0;
forAll(flux,energyJ)
{
	scatteringSource[energyI]+=sigmaS[energyI][energyJ][0] * flux[energyJ];
}

