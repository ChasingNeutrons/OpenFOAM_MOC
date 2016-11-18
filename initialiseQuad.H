/*
Create quadrant azimuthal angles, use angles to calculate number of tracks
originating from bottom and left boundary given spacing, geometry width 
and height.
This is then used to recalculate the effective number of angles thus giving the
effective track spacing
*/

//only calculate for quadrant - can apply angles and weights by symmetry
//Readjust azimuthal angles to make divisible by 2 and 4 - avoid cos(pi/2) etc.
if (naz%4!=0) {naz+= 4 - naz%4;}
label n2=naz/2;

List<scalar> phi(n2);		//List of azimuthal angles
List<label> nx(n2);		//List of numbers of tracks in x direction for angle i
List<label> ny(n2);		//List of numbers of tracks in y direction for angle i
List<scalar> tsi(n2);		//List of corrected track spacings for each angle i
label nnx, nny;
scalar p;

for(label i=0; i<naz/4; i++)
{
	label j= n2 - 1 - i;
	p=2*PI*(i+0.5)/naz;

	nnx=floor(w*mag(Foam::sin(p))/ts) + 1;
	nny=floor(h*mag(Foam::cos(p))/ts) + 1;

	nx[i]=nnx;
	nx[j]=nnx;
	ny[i]=nny;
	ny[j]=nny;
	
	
	//Calculate effective angle
	p=Foam::atan((h*nnx)/(w*nny));
	phi[i]=p;
	phi[j]=PI - p;
		
	tsi[i]=w*mag(Foam::sin(p))/nnx;
	tsi[j]=tsi[i];
}

//Create azimuthal quadrature
List<scalar> waz(n2);

forAll(phi, i)
{
	scalar w1, w2;
	if (i==0)
	{
		w1=phi[i];
		w2=(phi[i+1]-phi[i])/2;
	}
	else if (i==n2-1)
	{
		w1=(phi[i]-phi[i-1])/2;
		w2=PI - phi[i];
	}
	else
	{
		w1=(phi[i]-phi[i-1])/2;
		w2=(phi[i+1]-phi[i])/2;
	}
	waz[i]=w1+w2;
}

//Create weight and width quantity indexed by angle
List<scalar> weight_width(n2);
forAll(weight_width, i)
{
	weight_width[i]=waz[i]*tsi[i];
}

//Polar Quadrature - TY
List<scalar> wsintheta(npo);
List<scalar> sintheta(npo);
if (npo==1)
{
	wsintheta[0]= 2*0.798184;
	sintheta[0]=0.798184;
}else if (npo==3)
{
	wsintheta[0]=2.0*0.046233*0.166648;
        wsintheta[1]= 2.0*0.283619*0.537707;
        wsintheta[2]= 2.0*0.670148*0.932954;
	sintheta[0]=0.166648;
        sintheta[1]= 0.537707;
        sintheta[2]= 0.932954;
}else 
{
	wsintheta[0]=2.0*0.212854*0.3639000;
        wsintheta[1]= 2.0*0.787146*0.899900;
	sintheta[0]=0.3639000;
        sintheta[1]= 0.899900;
	if (npo != 2) {Info<< "polar angle must be an integer between 1 and 3 - polar angle has been reset to 2" << endl;}
}		

//Require costheta for anisotropic calculations
List<scalar>costheta(npo);
for(label j=0; j<npo; j++)
{
	costheta[j]=Foam::sqrt(1-sintheta[j]*sintheta[j]);
}
		
