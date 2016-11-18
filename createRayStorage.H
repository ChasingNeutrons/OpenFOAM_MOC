//Store ray info for each: [azimuthal angle] [ray in 4 cycle] [originating point along x or y axis]

//Make an array for storing BC?? Necessary??

//Store start and end point for each ray
int totalRays=0;
forAll(phi, i)
{
	totalRays+=nx[i]+ny[i];
}
List<List<point> > rayPoints(2);			//store ray start and end points
List<List<scalar> > alpha(2);				//albedo of boundary
List<int> raySegments(totalRays);			//Number of segments for each rays
List<DynamicList<scalar> > segLengths(totalRays);	//Store segment lengths for each ray
List<DynamicList<label> > cellIndices(totalRays);	//Store all cells traversed by a ray
List<scalar> rayAngle(totalRays);			//Store angle of ray - use for identifying complementary rays
List<label> angleInd(totalRays);			//Store ray angle index
List<scalar> wgta(totalRays);				//Store ray weight
forAll(rayPoints, i)
{
	rayPoints[i]=List<point>(totalRays);
	alpha[i]=List<scalar>(totalRays);
}

//For each angle, tally area of each cell
int totalCells=mesh.cells().size();
List<List<scalar> > approxArea(n2);	//Tally area approximations for each ray segment
forAll(approxArea, i)
{
	approxArea[i]=List<scalar>(totalCells);
	forAll(approxArea[i], c)
	{
		approxArea[i][c]=0.0;
	}
}

//For linear source problems store segment midpoints and entry points
List<List<DynamicList<scalar> > > xInCell(totalRays);	//Store all cell x-entrance points for a given ray
List<List<DynamicList<scalar> > > yInCell(totalRays);	//Store all cell y-entrance points for a given ray
forAll(xInCell, i)
{
	xInCell[i]=List<DynamicList<scalar> >(2);
	yInCell[i]=List<DynamicList<scalar> >(2);
}
List<DynamicList<scalar> > xCCell(totalRays);	//Store all cell x-centres for a given ray
List<DynamicList<scalar> > yCCell(totalRays);	//Store all cell y-centres for a given ray

//Store all sines and cosines for anisotropic calculations (multiply by -1 for opposite direction)
List<scalar>sinInd(n2);
List<scalar>cosInd(n2);
forAll(phi, i)
{
	sinInd[i]=Foam::sin(phi[i]);
	cosInd[i]=Foam::cos(phi[i]);
}

