//Transport Sweep Algorithm with a linear source
int sweeps=0;
do{
	//Reinitialise values
        err2=0.0;
        forAll(flux, energyI)
	{
		flux[energyI]*=0.0;
		fluxX[energyI]*=0.0;
		fluxY[energyI]*=0.0;
	}
		
	for(int lines=0; lines<klines; lines++)
	{	
		wa=wgta[lines];					//width and weight of line in azimuthal direction
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
		//Info<<sina<<endl;

		//List of areas approximated by lines at given angle
		List<scalar> approxAreaList=approxArea[angleInd[lines]];
		
		forAll(flux, energyI)
		{					
			for(int j=0; j<npo; j++)
			{
				wt = wa*wsintheta[j];	//combined polar weight and sintheta with azimuthal
				sinT=sintheta[j];
				List<scalar> f1 = F1[energyI][lines][j];
				List<scalar> hFac=H[energyI][lines][j];
	
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
						//Should these functions be *-1?
						sina*=-1;
						cosa*=-1;
					}

					//Take boundary flux for first cell
					afluxIn=aflux[energyI][j][lines][d];
					
					for(int k=c1; k!=c2+m; k+=m)
					{
						label ind=indList[k];
						A=area[ind];
						sigT=sigmaT[energyI][ind];
						scalar w_ASigT=wt/(A*sigT);
						scalar s=lengthList[k];
						scalar normFac=A/approxAreaList[ind];
						tau=sigT*s*normFac/sinT;

						//Confirmed that commenting fluxX,Y and relevant terms yields FS
						//Error in hFac, qh, qX, qY, fluxX, fluxY, Cx, Cy, Cxy, Mxx, Mxy, Myy				

						//Calculate expansion coefficients
						qb=one_over_4_PI*(Q[energyI][ind]+ xC[k]*qX[energyI][ind] 
							+ yC[k]*qY[energyI][ind]);

						qh=one_over_4_PI*(sinT*cosa*qX[energyI][ind] 
							+ sinT*sina*qY[energyI][ind])*normFac;
								
						//Calculate change in angular flux
						delta=(afluxIn - qb/sigT) * f1[k]
							-qh*(2*(tau-f1[k])-tau*f1[k])/(2*sigT*sigT);

						//Increment scalar flux and flux moments
						flux[energyI][ind] += delta*w_ASigT;
						fluxX[energyI][ind] += w_ASigT*(cosa*s*afluxIn*hFac[k]+xIn[k]*delta);
						fluxY[energyI][ind] += w_ASigT*(sina*s*afluxIn*hFac[k]+yIn[k]*delta);

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


	forAll(flux, energyI)
	{		
		//Increment flux moments by source moments		
		fluxX[energyI].internalField()+=(qX[energyI].internalField()*Cxx[energyI].internalField() 
			+ qY[energyI].internalField()*Cxy[energyI].internalField())/sigmaT[energyI].internalField();
		fluxY[energyI].internalField()+=(qY[energyI].internalField()*Cyy[energyI].internalField() 
			+ qX[energyI].internalField()*Cxy[energyI].internalField())/sigmaT[energyI].internalField();

		flux[energyI].internalField() += Q[energyI].internalField()/sigmaT[energyI].internalField();
		if(min(mag(flux[energyI])).value()<1e-6)
		{e1 = max(mag(flux[energyI]-prevFlux[energyI])).value();}
		else{e1 = max(mag(flux[energyI]-prevFlux[energyI])/flux[energyI]).value();}
		if( e1>err2 ) {err2=e1;}
		prevFlux[energyI]=flux[energyI];

		//Check error in angular flux
		forAll(aflux[energyI], j)
		{
			forAll(aflux[energyI][j], lines)
			{
				forAll(aflux[energyI][j][lines], d)
				{
					e1=mag(aflux[energyI][j][lines][d]-oldAflux[energyI][j][lines][d]);
					if(e1>err2) {err2=e1;}
					oldAflux[energyI][j][lines][d]=aflux[energyI][j][lines][d];
				}		
			}
		}
	}
	
        sweeps++;

}while( err2>1e-7);
Info<<"Transport sweep converged after "<<sweeps<<" iterations"<<endl;
