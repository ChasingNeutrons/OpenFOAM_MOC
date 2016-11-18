//Pre-compute and store Cgix, Cgiy and Cgixy

forAll(flux, energyI)
{
	for(int lines=0; lines<klines; lines++)
	{	
		wa=wgta[lines];
		List<scalar> xcList=xCCell[lines];
		List<scalar> ycList=yCCell[lines];		
		List<label> indList=cellIndices[lines];
		segNum=raySegments[lines];
		List<scalar> lengthList=segLengths[lines]; //list of azimuthal ray lengths
		label i0=angleInd[lines];			//index of azimuthal angle

		//Set trigonometric functions of azimuthal angle
		//Should these be strictly positive?
		sina=sinInd[i0];
		cosa=cosInd[i0];

		List<scalar> yIn;				//set list of local y-entry co-ords
		List<scalar> xIn;				//set list of local x-entry co-ords

		//List of areas approximated by lines at given angle
		List<scalar> approxAreaList=approxArea[angleInd[lines]];
		for(int j=0; j<npo; j++)
		{
			w = wa*wsintheta[j];
			sinT=sintheta[j];
			List<scalar> g1=G1[energyI][lines][j];
			List<scalar> g2=G2[energyI][lines][j];
			label i=i0;
	
			for(int d=0; d<2; d++) //are both directions necessary? a/2 index in reference!
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

				for(int k=c1; k!=c2+m; k+=m)
				{
					ind=indList[k];
					sigT=sigmaT[energyI][ind];
					scalar s=lengthList[k];
					A=area[ind];
					scalar t=s*A/approxAreaList[ind];
					scalar xc=xcList[k];
					scalar yc=ycList[k];

					Cm[energyI][i][j][ind]=g1[k]*t/sigT;						
					Cmx[energyI][i][j][ind]=(xc*g1[k] + cosa*s*g2[k]/2)*t/sigT;
					Cmy[energyI][i][j][ind]=(yc*g1[k] + sina*s*g2[k]/2)*t/sigT;

					//Calculate M factors if first energy iteration
					//integral of x^2 etc. along track with weighting factor
					// x= xin + cosa*t	y= yin + sina*t	where cosa and sina can be negative
					//int(x^2)=xin^2*t+cosa^2*t^3/3+xin*cosa*t^2
					//note division by two only for Mxy
					//remember division by sinT for track length!
					if(energyI==0)
					{
						//For calculating M
						scalar tm=t/sinT;
						scalar xi=xIn[k];
						scalar yi=yIn[k];
							
						Mxx[ind] += w*(xi*xi*tm+cosa*cosa*tm*tm*tm/3 + xi*cosa*tm*tm)/A;
						Myy[ind] += w*(yi*yi*tm+sina*sina*tm*tm*tm/3 + yi*sina*tm*tm)/A;
						Mxy[ind] += w*(xi*yi*tm+sina*cosa*tm*tm*tm/3 + (xi*sina+yi*cosa)*tm*tm/2)/A;
					}

				}
			}
		}
	}
}
