forAll(qm, energyI)
{
	forAll(qm[energyI], i)
	{
		if(i>=n2){i0=i-n2;}
		else{i0=i;}
		p=phi[i0];				//azimuthal angle required for spherical harmonic
		if(i>=n2){p+=PI;}

		forAll(qm[energyI][i], j)
		{
			qm[energyI][i][j] *= 0;
			mu=costheta[j];

			for(label energyJ=0; energyJ<energyGroups; energyJ++)
			{
				//Add all flux moments to fission or just zeroth moment?
				qm[energyI][i][j] += one_over_4_PI*chi[energyI] * nuSigmaEff[energyJ] * flux[energyJ]/keff;
				qm[energyI][i][j] += one_over_4_PI*sigmaS[energyI][energyJ][0]*flux[energyJ];
				for(label l=0; l<anisotropy+1; l++)
				{
					for(label r=-l; r<l+1; r++)
					{
						//Include factor of 2 to account for only covering half of polar angles
						//scalar R=2*sphericalHarmonic(l, r, p, mu);
						//qm[energyI][i][j] += R * sigmaS[energyI][energyJ][l] * fluxMo[energyJ][l][r+l];
						
					}
				}
				
			}
		}
	}
}
