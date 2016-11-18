//Found http://www.cfd-online.com/Forums/openfoam-programming-development/127449-ray-tracing-openfoam.html

//Set ray tracing options
treeBoundBox allBb(mesh.points());
   
scalar bbTol = 1e-12 * allBb.avgDim();
    
point& bbMin = allBb.min();
bbMin.x() -= bbTol;
bbMin.y() -= bbTol;
bbMin.z() -= bbTol;
    
point& bbMax = allBb.max();
bbMax.x() += 2*bbTol;
bbMax.y() += 2*bbTol;
bbMax.z() += 2*bbTol;
    
indexedOctree<treeDataFace> faceTree
(
	treeDataFace(false, mesh),       //formerly mesh_
	allBb, // overall search domain
	8, // maxLevel
	10, // leafsize
	3.0 // duplicity
);

label klines=0;	//track number of rays - index into storage for each ray

//Loop over all azimuthal angles
forAll(phi, i)							
{
	scalar t=tsi[i];	
	scalar p=phi[i];
	point pStart;
	point pBegin;
	point pEnd;
	int inc;
	scalar dx=t/mag(Foam::sin(p));
	scalar dy=t/mag(Foam::cos(p));
	vector dxvec(dx, 0.0, 0.0);
	vector dyvec(0.0, dy, 0.0);

	if(phi[i]<PI/2.0)
	{
		pStart=vector(0.0, 0.0, d_over_2);	//BOTTOM LEFT 	//point from which to begin all ray tracing
		inc=+1;					//direction in which rays are incremented on the x-axis
	}else
	{
		pStart=vector(w, 0.0, d_over_2);	//BOTTOM RIGHT 	//point from which to begin all ray tracing
		inc=-1;					//direction in which rays are incremented on the x-axis
	}
	
	label nnx=nx[i];
	label nny=ny[i];
	label nmax=nnx + nny;
	DynamicList<label> angleLine;			//list to store lines corresponding to a given angular index

	for (int n=0; n<nmax; n++)
	{		
		//Begin ray tracing from each starting point on x and y axis
		if(n<nnx)
		{
			pBegin=pStart + inc*(n+0.5)*dxvec; 	//plot along x-axis starting from origin +0.5*x-spacing
		}else
		{
			pBegin=pStart + (n-nnx+0.5)*dyvec;	//plot along y-axis starting from origin +0.5*y-spacing
		}
		
		scalar x0;
		if (p<PI/2.0){x0=w;}
		else{x0=0.0;}
		
		scalar ytest=Foam::tan(p)*(x0-pBegin.x()) + pBegin.y();
		
		if(ytest<h)
		{
			pEnd=vector(x0, ytest, d_over_2);
		}else
		{
			pEnd=vector( ((h-pBegin.y())/Foam::tan(p) + pBegin.x()), h, d_over_2);
		}
		
		//For ray klines, set start and end point
		rayPoints[0][klines]=pBegin;
		rayPoints[1][klines]=pEnd;
		
		//Set boundary conditions based on start and endpoint
		for (int r=0; r<2; r++)
		{
			point pRay=rayPoints[r][klines];
			scalar yval=pRay.y();
			scalar xval=pRay.x();
			//Info<<pRay<<endl;
			
			//Direction of ray hitting boundary
			int db;
			//scalar p2;
			//if(p<PI/2){p2 = p;}
			//else {p2 = PI-p;}

			if(r==0){db=1;}
			else{db=0;}
			if(mag(xval)<1e-12) 
			{
				alpha[db][klines]=alphaL;
				//leakFactor[db][klines]=t*waz[i]*Foam::cos(p2);
			}
			else if(mag(xval-w)<1e-12) 
			{
				alpha[db][klines]=alphaR;
				//leakFactor[db][klines]=t*waz[i]*Foam::cos(p2);
			}
			else if(mag(yval)<1e-12) 
			{
				alpha[db][klines]=alphaB;
				//leakFactor[db][klines]=t*waz[i]*Foam::sin(p2);
			}
			else if(mag(yval-h)<1e-12) 
			{
				alpha[db][klines]=alphaT;
				//leakFactor[db][klines]=t*waz[i]*Foam::sin(p2);
			}
			else {Info<<"Error in identifying boundary condition"<<endl;}
		}
		
		//Store ray angle
		rayAngle[klines]=p;
	   
		vector eVec(pEnd - pBegin); // line vector
			    
		//const scalar eMag = mag(eVec); // edge length
		vector tolVec = 1e-12*eVec;    
		point p0 = pBegin + tolVec;
		point p1 = pEnd + tolVec;

		//Storage for all segment lengths and cell labels
		DynamicList<scalar> Lray;
		DynamicList<int> cellLabel;
		DynamicList<scalar> rayArea;
		DynamicList<scalar> xInList;
		DynamicList<scalar> xOutList;
		DynamicList<scalar> yInList;
		DynamicList<scalar> yOutList;
		DynamicList<scalar> xCList;
		DynamicList<scalar> yCList;
		int k=0;

		while(true)
		{
			pointIndexHit pHit = faceTree.findLine(p0, p1);
			if (pHit.hit())
			{
				//label faceI = pHit.index(); // face label of hit face - can be used to identify boundary!
				point hitPoint = pHit.hitPoint(); // intersection point

				scalar length=mag(hitPoint-p0);
				label c=mesh.findCell(p0);
				scalar xCentre=centres[c].component(0);
				scalar yCentre=centres[c].component(1);

				cellLabel.append(c);			//store cell indices passed through
				Lray.append(length);			//store ray lengths
				approxArea[i][c]+=length*t;		//increment approximate cell area
				
				//For linear source problems store local co-ordinates of cell entry point and ray centres
				xInList.append(p0.x()-xCentre);
				xOutList.append(hitPoint.x()-xCentre);
				xCList.append((p0.x()+hitPoint.x())/2 - xCentre);
				yInList.append(p0.y()-yCentre);
				yOutList.append(hitPoint.y()-yCentre);
				yCList.append((p0.y()+hitPoint.y())/2 - yCentre);
				
				vector area = mesh.faceAreas()[pHit.index()];
				scalar typDim = Foam::sqrt(mag(area));
			
				//Increment number of cells searched
				k++;
					
				// stop search if new start point is near to edge end
				if ((mag(hitPoint - pEnd)/typDim) < 1e-6)
				{
					//Info<<Lray<<endl;					
					raySegments[klines]=k;			//store number of ray segments
					segLengths[klines].append(Lray);	//store length of segments
					cellIndices[klines].append(cellLabel);	//store cell labels
					wgta[klines]=waz[i]*t;			//store ray weight
					rayAngle[klines]=p;			//store ray angle
					angleInd[klines]=i;			//store index of ray angle
					xInCell[klines][0].append(xInList);
					xInCell[klines][1].append(xOutList);
					yInCell[klines][0].append(yInList);
					yInCell[klines][1].append(yOutList);
					xCCell[klines].append(xCList);
					yCCell[klines].append(yCList);
					klines++;
					break;
				}

				// set new start point shortly after previous start point
				p0 = hitPoint + tolVec;
				label newCell=mesh.findCell(p0);
				//Ensure p0 is not in the same cell!
				while(newCell==c){
					p0+=tolVec;
					newCell=mesh.findCell(p0);
				} 
			}
			else
			{	
				// No hit.
				scalar length=mag(p1-p0);
				label c=mesh.findCell(p0);
				scalar xCentre=centres[c].component(0);
				scalar yCentre=centres[c].component(1);

				cellLabel.append(c);			//store cell indices passed through
				Lray.append(length);			//store ray lengths
				approxArea[i][c]+=length*t;		//increment approximate cell area

				//For linear source problems store local co-ordinates of cell entry point
				xInList.append(p0.x()-xCentre);
				xOutList.append(p1.x()-xCentre);
				xCList.append((p0.x()+p1.x())/2 - xCentre);
				yInList.append(p0.y()-yCentre);
				yOutList.append(p1.y()-yCentre);
				yCList.append((p0.y()+p1.y())/2 - xCentre);
				
				raySegments[klines]=k;			//store number of ray segments
				segLengths[klines].append(Lray);	//store length of segments
				cellIndices[klines].append(cellLabel);	//store cell labels
				wgta[klines]=waz[i]*t;			//store ray weight
				rayAngle[klines]=p;			//store ray angle
				xInCell[klines][0].append(xInList);
				xInCell[klines][1].append(xOutList);
				yInCell[klines][0].append(yInList);
				yInCell[klines][1].append(yOutList);
				xCCell[klines].append(xCList);
				yCCell[klines].append(yCList);
				klines++;				
				Info<<"Ray tracing failure"<<endl;
				Info<<"Line "<<klines-1<<" Start point = "<<pBegin<<endl;
				Info<<"Angle ="<<p<<endl;
				Info<<p0<<endl;
				Info<<pEnd<<endl;
				break;
			}
		}
	}
}
