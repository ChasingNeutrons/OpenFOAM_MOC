//Legendre polynomial function up to fifth order
scalar legendrePolynomial(label l, label r, scalar mu)
{
	if(l==0){ return 1.0;}
	else if(l==1)
	{
		if(r==0){return mu;}
		else if(r==1){ return -Foam::sqrt(1.0-mu*mu);}
		else{
			Info<<"Incorrect Legendre input"<<endl;
			return EXIT_FAILURE;
		}

	}else if(l==2)
	{
		if(r==0){return 0.5*(3.0*mu*mu-1.0);}
		else if(r==1){return -3.0*mu*Foam::sqrt(1.0-mu*mu);}
		else if(r==2){return 3.0*(1.0-mu*mu);;}
		else{
			Info<<"Incorrect Legendre input"<<endl;
			return EXIT_FAILURE;
		}
		
	}else if(l==3)
	{
		if(r==0){return 0.5*(5.0*mu*mu*mu - 3.0*mu);}
		else if(r==1){return -1.5*(5*mu*mu-1.0)*Foam::sqrt(1.0-mu*mu);}
		else if(r==2){return 15.0*mu*(1.0-mu*mu);}
		else if(r==3) {return -15.0*Foam::pow((1.0-mu*mu),1.5);}
		else{
			Info<<"Incorrect Legendre input"<<endl;
			return EXIT_FAILURE;
		}

	}else if(l==4)
	{
		if (r==0){return (35.0*mu*mu*mu*mu-30*mu*mu+3.0)/8.0;}
		else if(r==1){return -2.5*(7.0*mu*mu*mu-3*mu)*Foam::sqrt(1.0-mu*mu);}
		else if(r==2){return 7.5*(7*mu*mu-1)*(1-mu*mu);}
		else if(r==3){return -105.0*mu*Foam::pow((1-mu*mu),1.5);}
		else if(r==4){return 105.0*Foam::sqr(1-mu*mu);}
		else{
			Info<<"Incorrect Legendre input"<<endl;
			return EXIT_FAILURE;
		}
	}else{
		Info<<"Too high a Legendre order has been entered"<<endl;		
		return EXIT_FAILURE;
	}

}
	


//Function to evaluate spherical harmonic
scalar sphericalHarmonic(label l, label r, scalar phi, scalar mu)
{
	if(l==0){return Foam::sqrt(one_over_4_PI);}
	scalar normalisation=Foam::sqrt((2*l+1)*one_over_4_PI);
	scalar fac;
	scalar trig;
	label magr=mag(r);
	if(r==0){fac=1.0; trig=1.0;}
	else
	{
		fac=Foam::sqrt(2.0*Foam::factorial(l-magr)/Foam::factorial(l+magr));
		if(r>0){trig=Foam::cos(r*phi);}
		else{trig=Foam::sin(magr*phi);}
	}

	return normalisation*fac*legendrePolynomial(l, magr, mu)*trig;
}
