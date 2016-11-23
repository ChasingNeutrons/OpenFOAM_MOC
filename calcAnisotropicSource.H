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
				//Fission source = one_over_pi*nuSigmaEff*integral(psi)*d(omega)
				//integral(psi)*d(omega) = normal scalar flux
				qm[energyI][i][j] += one_over_4_PI * chi[energyI] * nuSigmaEff[energyJ] * flux[energyJ]/keff;
				
				for(label l=0; l<anisotropy+1; l++)
				{
					for(label r=-l; r<l+1; r++)
					{
						//Scattering source = Rlr(phi,theta)*sigmaSl*fluxlr
						scalar R1=sphericalHarmonic(l, r, p, mu);
						scalar R2=sphericalHarmonic(l, r, p, -mu);
						qm[energyI][i][j] += (R1+R2) * sigmaS[energyI][energyJ][l] * fluxMo[energyJ][l][r+l]/2.0;
						
					}
				}				
			}
		}
	}
}
