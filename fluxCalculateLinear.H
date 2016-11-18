//Reinitialise flux
forAll(fluxMo, energyI)
{
	for(label l=0; l<anisotropy+1; l++)
	{
		for(label r=-l; r<l+1; r++)
		{
			fluxMo[energyI][l][r+l] *=0.0;
			fluxXMo[energyI][l][r+l] *=0.0;
			fluxYMo[energyI][l][r+l] *=0.0;
		}
	}
}

//Total up scalar fluxes from previously computed values of average angular fluxes
label mult=1;			//constant for modifying sin and cos based on azimuthal direction
for(label i=0; i<naz; i++)
{
	if(i>=n2)
	{
		i0=i-n2;
		mult=-1;
	}
	else{
		i0=i;
		mult=1;
	}	
	wa=weight_width[i0];				//width and weight of line in azimuthal direction
	sina=mult*sinInd[i0];
	cosa=mult*cosInd[i0];
	p=phi[i0];				//azimuthal angle required for spherical harmonic
	if(i>=n2){p+=PI;}
	scalarField normField=area/approxArea[i0];	//area normalisation for volume conservation
		
	forAll(fluxMo, energyI)
	{					
		for(label j=0; j<npo; j++)
		{
			wt = wa*wsintheta[j];	//combined polar weight and sintheta with azimuthal
			sinT=sintheta[j];	
			mu=costheta[j];		//polar cosine for calculating spherical harmonic

			//Increment scalar flux and flux moments
			//Can sigT cancel out? Or include sigT in summation ???
			for(label l=0; l<anisotropy+1; l++)
			{
				for(label r=-l; r<l+1; r++)
				{
					scalar R=sphericalHarmonic(l,r,p,mu);

					fluxMo[energyI][l][r+l].internalField()+=psi[energyI][i][j].internalField()*wt*R/area;

					fluxXMo[energyI][l][r+l].internalField()+=wt*R*(sinT*cosa*psiH[energyI][i][j].internalField()/normField
							+ psiX[energyI][i][j].internalField())/area;

					fluxYMo[energyI][l][r+l].internalField()+=wt*R*(sinT*sina*psiH[energyI][i][j].internalField()/normField
							+ psiY[energyI][i][j].internalField())/area;
				}
			}
		}
	}	
}

