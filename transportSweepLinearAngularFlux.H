//Transport Sweep Algorithm with a linear source
int sweeps=0;
do{
	//Reinitialise average angular flux
        err2=0.0;
        forAll(psi, energyI)
	{
		forAll(psi[energyI], i)
		{
			forAll(psi[energyI][i], j)
			{

				psi[energyI][i][j] *=0.0;
				psiX[energyI][i][j]*=0.0;
				psiY[energyI][i][j]*=0.0;
				psiH[energyI][i][j]*=0.0;
			}
		}
	}
		
	for(int lines=0; lines<klines; lines++)
	{	
		List<label> indList=cellIndices[lines];		//list of indices of cells intersected
		segNum=raySegments[lines];			//number of segments of ray
		List<scalar> yIn;				//set list of local y-entry co-ords
		List<scalar> xIn;				//set list of local x-entry co-ords
		List<scalar> yC=yCCell[lines];
		List<scalar> xC=xCCell[lines];
		List<scalar> lengthList=segLengths[lines];	//list of azimuthal ray lengths
		label i0=angleInd[lines];			//index of azimuthal angle

		//Set trigonometric functions of azimuthal angle
		//Should these be strictly positive? (if so make sure to remove *-1 in the direction loop)
		sina=sinInd[i0];
		cosa=cosInd[i0];

		//List of areas approximated by lines at given angle
		List<scalar> approxAreaList=approxArea[angleInd[lines]];
		
		forAll(flux, energyI)
		{					
			for(int j=0; j<npo; j++)
			{
				sinT=sintheta[j];
				List<scalar> f1 = F1[energyI][lines][j];
				List<scalar> hFac=H[energyI][lines][j];
				label i=i0;
	
				for(int d=0; d<2; d++)
				{
					yIn=yInCell[lines][d];
					xIn=xInCell[lines][d];
					if (d==0)
					{
						c1=0;
						c2=segNum-1;
						m=+1;
					}else if (d==1)
					{
						c1=segNum-1;
						c2=0;
						m=-1;
						sina*=-1;
						cosa*=-1;
						i=i0+n2;
					}

					//Take boundary flux for first cell
					afluxIn=aflux[energyI][j][lines][d];
							
					for(int k=c1; k!=c2+m; k+=m)
					{
						label ind=indList[k];
						A=area[ind];
						sigT=sigmaT[energyI][ind];
						scalar s=lengthList[k];
						scalar normFac=A/approxAreaList[ind];
						tau=sigT*s*normFac/sinT;
						scalar t=s*normFac;
						scalar xi=xIn[k];
						scalar yi=yIn[k];
						scalar one_x_y=1+xi+yi;
						scalar q=qm[energyI][i][j][ind];
						scalar qx=qmX[energyI][i][j][ind];
						scalar qy=qmY[energyI][i][j][ind];

						//Calculate expansion coefficients
						//Should the 4 pi still be here for the angular source?
						qb=one_over_4_PI*(q + xC[k]*qx + yC[k]*qy);

						qh=one_over_4_PI*(sinT*cosa*qx	+ sinT*sina*qy)/normFac;
								
						//Calculate change in angular flux
						delta=(afluxIn - qb/sigT) * f1[k]
							-qh*(2*tau-(2+tau)*f1[k])/(2*sigT*sigT);

						//Increment average angular flux and flux moments
						psi[energyI][i][j][ind] += delta/sigT 
							+ t*one_x_y*q*one_over_4_PI/(sinT*sigT);

						psiX[energyI][i][j][ind] += xi*delta/sigT
							+ t*xC[k]*one_x_y*q*one_over_4_PI/(sinT*sigT);

						psiY[energyI][i][j][ind] += yi*delta/sigT
							+ t*yC[k]*one_x_y*q*one_over_4_PI/(sinT*sigT);

						//Increment psi hat
						psiH[energyI][i][j][ind]+=afluxIn*hFac[k]*t/(sinT*sigT);

						//Calculate outgoing angular flux
						afluxIn -= delta;

						//Apply boundary conditions if exiting boundary cell
						if (k==c2)
						{
							albedo=alpha[d][lines];				
							refRay=compRay[lines][d];
							refDir=compDir[lines][d];
							aflux[energyI][j][refRay][refDir]=albedo*afluxIn;
						}
					}
				}
			}
		}
	}

	//Increment psi hat by source terms
	forAll(psiH, energyI)
	{
		forAll(psiH[energyI], i)
		{
			forAll(psiH[energyI][i], j)
			{
				psiH[energyI][i][j].internalField() += 
					one_over_4_PI*(qm[energyI][i][j].internalField() * Cm[energyI][i][j].internalField() +
					qmX[energyI][i][j].internalField() * Cmx[energyI][i][j].internalField() +
					qmY[energyI][i][j].internalField() * Cmy[energyI][i][j].internalField())
					/(sintheta[j]*sigmaT[energyI].internalField());
			}
		}
	}

	//Check error
	forAll(aflux, energyI)
	{		
		forAll(aflux[energyI], j)
		{
			forAll(aflux[energyI][j], lines)
			{
				forAll(aflux[energyI][j][lines], d)
				{
						e1=mag(aflux[energyI][j][lines][d]-oldAflux[energyI][j][lines][d])/aflux[energyI][j][lines][d];
						if(e1>err2) {err2=e1;}
						oldAflux[energyI][j][lines][d]=aflux[energyI][j][lines][d];
				}
			}			
		}
	}		
        sweeps++;

}while( err2>1e-7);
Info<<"Transport sweep converged after "<<sweeps<<" iterations"<<endl;
