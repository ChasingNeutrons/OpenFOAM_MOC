List<List<List<List<scalar> > > > F1(energyGroups);
   
forAll(F1, energyI)
{
	F1[energyI]=List<List<List<scalar> > >(klines);
	
	forAll(F1[energyI], lines)
	{
		F1[energyI][lines]=
		List<List<scalar> >(npo);
		
		
		forAll(F1[energyI][lines], j)
		{
			F1[energyI][lines][j]=
			List<scalar>(raySegments[lines]);

			forAll(F1[energyI][lines][j], seg)
			{
				label ind = cellIndices[lines][seg];				
				
				F1[energyI][lines][j][seg]= 
				1.0 -Foam::exp(-segLengths[lines][seg] 
				* sigmaT[energyI][ind]/sintheta[j]);
			}
		}
	}
}

List<List<List<List<scalar> > > > H(energyGroups);
List<List<List<List<scalar> > > > G1(energyGroups);
List<List<List<List<scalar> > > > G2(energyGroups);
if(linearSource)
{	   
	forAll(H, energyI)
	{
		H[energyI]=List<List<List<scalar> > >(klines);
		G1[energyI]=List<List<List<scalar> > >(klines);		
		G2[energyI]=List<List<List<scalar> > >(klines);
	
		forAll(H[energyI], lines)
		{
			H[energyI][lines]=
			List<List<scalar> >(npo);

			G1[energyI][lines]=
			List<List<scalar> >(npo);

			G2[energyI][lines]=
			List<List<scalar> >(npo);
		
		
			forAll(H[energyI][lines], j)
			{
				H[energyI][lines][j]=
				List<scalar>(raySegments[lines]);

				G1[energyI][lines][j]=
				List<scalar>(raySegments[lines]);				
				
				G2[energyI][lines][j]=
				List<scalar>(raySegments[lines]);

				forAll(H[energyI][lines][j], seg)
				{
					label ind = cellIndices[lines][seg];				
				
					scalar tau=segLengths[lines][seg] * sigmaT[energyI][ind]/sintheta[j];
					scalar f1=F1[energyI][lines][j][seg];
					scalar G=1 + tau/2 - (1 + 1/tau)*f1;
				
					H[energyI][lines][j][seg]=tau/2 - G;
					G1[energyI][lines][j][seg]=G;
					G2[energyI][lines][j][seg]=2*tau/3 -(1 + 2/tau)*G;
				}
			}
		}
	}
}
