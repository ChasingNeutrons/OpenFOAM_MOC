//Legendre polynomial function up to fifth order
scalar legendrePolynomial(label l, label r, scalar mu)
{
	scalar p0, p1, p2, c;
	label magr=mag(r);
	if(r<0){c=Foam::pow(-1,magr)*Foam::factorial(l-magr)/Foam::factorial(l+magr);}
	else{c=1;}

	label count=0;
	p1=1;
	p0=0;

	if(l==0){return p1;}
	
	//Find r=0 polynomial
	do{
		if(count==0){p2=mu;}
		else{
			p2=((2*count+1)*mu*p1-count*p0)/(count+1);
		}
		p0=p1;
		p1=p2;

		count++;
	}while(count<l);

	count=0;

	if(r==0){return p1;}

	//Iterate until find polynomial with correct order in r
	do{
		p2=((l-count)*mu*p1 - (l+count)*p0)/Foam::sqrt(1-mu*mu);
		p0=((l-count+1)*p1 - (l+count+1)*mu*p0)/Foam::sqrt(1-mu*mu);
		p1=p2;
		count++;
	}while(count<magr);

	return c*p1;

}
	


//Function to evaluate spherical harmonic
scalar sphericalHarmonic(label l, label r, scalar phi, scalar mu)
{
	scalar normalisation=Foam::sqrt((2*l+1)*one_over_4_PI);
	scalar fac;
	scalar trig;
	if(r==0){fac=1; trig=1;}
	else
	{
		fac=Foam::sqrt(2.0*Foam::factorial(l-r)/Foam::factorial(l+r));
		if(r>0){trig=Foam::cos(r*phi);}
		else{trig=Foam::sin(r*phi);}
	}

	return normalisation*fac*legendrePolynomial(l, r, mu)*trig;
}
