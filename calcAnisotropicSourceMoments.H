forAll(QmX, energyI)
{
	for(label i=0; i<naz; i++)
	{
		if(i>=n2){i0=i-n2;}
		else{i0=i;}
		p=phi[i0];				//azimuthal angle required for spherical harmonic
		if(i>=n2){p+=PI;}
		
		for(label j=0; j<npo; j++)
		{
	
			QmX[energyI][i][j] *= 0;
			QmY[energyI][i][j] *= 0;
			mu=costheta[j];		//polar cosine for calculating spherical harmonic
	
			//Add fission terms
			//Note, formula givn in paper does not include 4PI
			forAll(flux, energyJ)
			{
				QmX[energyI][i][j] += chi[energyI] * nuSigmaEff[energyJ] * fluxX[energyJ];
				QmY[energyI][i][j] += chi[energyI] * nuSigmaEff[energyJ] * fluxY[energyJ];
			
				//Include anisotropic scattering contributions
				for(label l=0; l<anisotropy+1; l++)
				{
					for(label r=-l; r<l+1; r++)
					{
						scalar R=sphericalHarmonic(l, r, p, mu);
						QmX[energyI][i][j] += R * sigmaS[energyI][energyJ][l] * fluxXMo[energyJ][l][r+l];
						QmY[energyI][i][j] += R * sigmaS[energyI][energyJ][l] * fluxYMo[energyJ][l][r+l];
					}
				}
			}
			
			//Calculate x and y expansion coefficients by inverting the governing linear system:
			//M*qxy = Qxy
			qmX[energyI][i][j].internalField() = (Myy.internalField()*QmX[energyI][i][j].internalField()-Mxy.internalField()*QmY[energyI][i][j].internalField())
						/(Myy.internalField()*Mxx.internalField() -Mxy.internalField()*Mxy.internalField());
			qmY[energyI][i][j].internalField() = (Mxx.internalField()*QmY[energyI][i][j].internalField()-Mxy.internalField()*QmX[energyI][i][j].internalField())
						/(Myy.internalField()*Mxx.internalField() -Mxy.internalField()*Mxy.internalField());
		}
	}
}

