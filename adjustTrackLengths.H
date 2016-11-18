//Readjust track lengths to account for area approximation
forAll(approxArea, lines)
{
	label azInd=angleInd[lines];	
	for(int k=0; k<raySegments[lines]; k++)
	{		
		label cellInd=cellIndices[lines][k];
		scalar cellArea=approxArea[azInd][cellInd];
		if(cellArea<=0)
		{
			Info<<"Track spacing may not be sufficiently tight"<<endl;
		}
		else
		{
			segLengths[lines][k]*=area[cellInd]/cellArea;
		}
	}
}