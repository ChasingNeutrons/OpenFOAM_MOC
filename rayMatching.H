//Find complementary tracks corresponding to ray travelling in both directions
List<List<label> > compRay(totalRays);	//label of complementary ray
List<List<label> > compDir(totalRays);	//orientation of complementary ray

forAll(compRay, lines)
{
	compRay[lines]= List<label>(2);
	compDir[lines]= List<label>(2);
}

scalar angle;
scalar angle2;
scalar compAngle;
point pRay1End;
point pRay2Start;
bool found;

for(int lines=0; lines<totalRays; lines++)
{
	angle=rayAngle[lines];
	compAngle=PI - angle;
	
	for(int dir=0; dir<2; dir++)
	{
		found=0;
		pRay1End=rayPoints[dir][lines];
	
		forAll(compRay, compLines)
		{
			angle2=rayAngle[compLines];
			
			if(mag(compAngle-angle2)<1e-12)
			{
				for(int dd=0; dd<2; dd++)
				{
					pRay2Start=rayPoints[dd][compLines];
					
					if(mag(pRay2Start-pRay1End)<1e-12)
					{
						compRay[lines][dir]=compLines;
						compDir[lines][dir]=dd;
						
						found=1;
						break;
					}
				}
			}if(found){break;}
		}
	}
}



				
