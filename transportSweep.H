//Transport Sweep Algorithm
int sweeps=0;
do{
	//Reinitialise values
        err2=0.0;
        forAll(flux, energyI)
	{
		flux[energyI]*=0.0;
	}
		
	for(int lines=0; lines<klines; lines++)
	{	
		wa=wgta[lines];
		List<label> indList=cellIndices[lines];
		segNum=raySegments[lines];
		List<scalar> lengthList=segLengths[lines];	//list of azimuthal ray lengths

		//List of areas approximated by lines at given angle
		List<scalar> approxAreaList=approxArea[angleInd[lines]];
		
		forAll(flux, energyI)
		{					
			for(int j=0; j<npo; j++)
			{
				wt = wa*wsintheta[j];
				sinT=sintheta[j];
				List<scalar> f1 = F1[energyI][lines][j];
	
				for(int d=0; d<2; d++)
				{

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
					}
					

					//Take boundary flux for first cell along ray
					afluxIn=aflux[energyI][j][lines][d];
							
					for(int k=c1; k!=c2+m; k+=m)
					{
						ind=indList[k];
						sigT=sigmaT[energyI][ind];
						A=area[ind];
						tau=sigT*lengthList[k]*A/approxAreaList[ind]/sinT;
								
						//Calculate change in angular flux
						delta=(afluxIn - one_over_4_PI*Q[energyI][ind]/sigT) * f1[k];

						//Increment scalar flux
						flux[energyI][ind] += delta*wt/(A*sigT);

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
