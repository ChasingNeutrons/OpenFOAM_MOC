//Reinitialise flux
forAll(fluxMo, energyI)
{
	flux[energyI] *= 0.0;
	for(label l=0; l<anisotropy+1; l++)
	{
		for(label r=-l; r<l+1; r++)
		{
			fluxMo[energyI][l][r+l] *=0.0;
		}
	}
}

//Total up scalar fluxes from previously computed values of average angular fluxes
scalarField fluxInc;
for(label i=0; i<naz; i++)
{
	if(i>=n2){i0=i-n2;}
	else{i0=i;}
	wa=weight_width[i0];			//width and weight of line in azimuthal direction
	p=phi[i0];				//azimuthal angle required for spherical harmonic
	if(i>=n2){p+=PI;}
		
	forAll(flux, energyI)
	{					
		for(label j=0; j<npo; j++)
		{
			wt = wa*wsintheta[j];	//combined polar weight and sintheta with azimuthal
			mu=costheta[j];		//polar cosine for calculating spherical harmonic

			fluxInc=psi[energyI][i][j].internalField()*wt/area;
			flux[energyI].internalField()+=fluxInc;

			//Increment scalar flux and flux moments
			for(label l=0; l<anisotropy+1; l++)
			{
				for(label r=-l; r<l+1; r++)
				{
					//Spherical harmonic used here is normalised to 1 i.e. int(R*R)d(omega)=1
					// philr = integral(psi * Ylr)d(omega)
					fluxMo[energyI][l][r+l].internalField()+=fluxInc*(sphericalHarmonic(l,r,p,mu)+sphericalHarmonic(l,r,p,-mu))/2.0;
				}
			}
		}
	}	
}
