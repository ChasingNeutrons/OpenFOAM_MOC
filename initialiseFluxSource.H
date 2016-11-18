List<List<List<List<scalar> > > > aflux(energyGroups);
List<List<List<List<scalar> > > > oldAflux(energyGroups);
   
Info<<"Initialising fluxes and sources"<<endl;
if(anisotropy>0)
{
	#include "calcAnisotropicSource.H"
}
forAll(flux, energyI)
{
	if(anisotropy==0)
	{
		#include "calcScatteringSource.H"
		#include "calcFissionSource.H"
		Q.set(energyI,(scatteringSource[energyI] + fissionSource[energyI]/keff));
		oldQ.set(energyI, Q[energyI]);
	}
	
	prevFlux.set(energyI, flux[energyI]);
	oldFlux.set(energyI, flux[energyI]);	
	
	aflux[energyI]=List<List<List<scalar> > >(npo);
	
	oldAflux[energyI]=List<List<List<scalar> > >(npo);
	
	forAll(aflux[energyI], j)
	{
		aflux[energyI][j]=List<List<scalar> >(totalRays);
	
		oldAflux[energyI][j]=List<List<scalar> >(totalRays);

		forAll(aflux[energyI][j], lines)
		{
			aflux[energyI][j][lines]=List<scalar>(2);

			oldAflux[energyI][j][lines]=List<scalar>(2);
			
			forAll(aflux[energyI][j][lines], d)
			{
				aflux[energyI][j][lines][d]=1.0;
				oldAflux[energyI][j][lines][d]=1.0;
			}
		}
	}
}

#include "calcTotalFission.H"
prevFissions=gSum(totalFissions.internalField() * area);
