
Info<< "Creating scalar flux\n" << endl;
volScalarField defaultFlux
(
	IOobject
	(
		"defaultFlux",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

//Create fluxes
PtrList<volScalarField> flux(energyGroups);
forAll(flux,fluxI)
{
       Info<< "    Adding to flux\n" << endl;
       flux.set(
		fluxI,
       		new volScalarField
		(		       
			IOobject
			(
			   "flux"+Foam::name(fluxI),
			    runTime.timeName(),
			    mesh,
			    IOobject::READ_IF_PRESENT,
			    IOobject::AUTO_WRITE
			),
			defaultFlux
		)
	);
}

PtrList<PtrList<PtrList<volScalarField> > >fluxMo(energyGroups);
forAll(fluxMo,fluxI)
{
       fluxMo.set(fluxI,new PtrList<PtrList<volScalarField> >(anisotropy+1));
       forAll(fluxMo[fluxI], l)
       {
		fluxMo[fluxI].set(l, new PtrList<volScalarField>(2*l+1));
		forAll(fluxMo[fluxI][l], r)
		{
			fluxMo[fluxI][l].set
			(
				r,
				new volScalarField
				(		       
					IOobject
					(
					   "fluxMo"+Foam::name(fluxI)+Foam::name(l)+Foam::name(r),
					    runTime.timeName(),
					    mesh,
					    IOobject::READ_IF_PRESENT,
					    IOobject::NO_WRITE
					),
					defaultFlux
				)
			);
		}
        }
}

//Create flux x moments
Info<< "Creating scalar flux x-moment\n" << endl;
PtrList<volScalarField> fluxX(energyGroups);

forAll(fluxX,fluxXI)
{
       Info<< "    Adding to flux x-moment\n" << endl;
       fluxX.set(
		fluxXI,
		new volScalarField
		(		       
			IOobject
			(
			   "fluxX"+Foam::name(fluxXI),
			    runTime.timeName(),
			    mesh,
			    IOobject::READ_IF_PRESENT,
			    IOobject::NO_WRITE
			),
		mesh,
		dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
		zeroGradientFvPatchScalarField::typeName
		)
	);
}

PtrList<PtrList<PtrList<volScalarField> > >fluxXMo(energyGroups);

forAll(fluxXMo,fluxXI)
{
       fluxXMo.set(fluxXI,new PtrList<PtrList<volScalarField> >(anisotropy+1));
       forAll(fluxXMo[fluxXI], l)
       {
		fluxXMo[fluxXI].set(l, new PtrList<volScalarField>(2*l+1));
		forAll(fluxXMo[fluxXI][l], r)
		{
			fluxXMo[fluxXI][l].set
			(
				r,
				new volScalarField
				(		       
					IOobject
					(
					   "fluxXMo"+Foam::name(fluxXI)+Foam::name(l)+Foam::name(r),
					    runTime.timeName(),
					    mesh,
					    IOobject::READ_IF_PRESENT,
					    IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
        }
}

//Create flux y moments
Info<< "Creating scalar flux y-moment\n" << endl;
PtrList<volScalarField> fluxY(energyGroups);

forAll(fluxY,fluxYI)
{
       Info<< "    Adding to flux y-moment\n" << endl;
       fluxY.set(fluxYI,
		new volScalarField
		(		       
			IOobject
			(
			   "fluxY"+Foam::name(fluxYI),
			    runTime.timeName(),
			    mesh,
			    IOobject::READ_IF_PRESENT,
			    IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

PtrList<PtrList<PtrList<volScalarField> > >fluxYMo(energyGroups);

forAll(fluxYMo,fluxYI)
{
       Info<< "    Adding to flux y-moment\n" << endl;
       fluxYMo.set(fluxYI,new PtrList<PtrList<volScalarField> >(anisotropy+1));
       forAll(fluxYMo[fluxYI], l)
       {
		fluxYMo[fluxYI].set(l, new PtrList<volScalarField>(2*l+1));
		forAll(fluxYMo[fluxYI][l], r)
		{
			fluxYMo[fluxYI][l].set
			(
				r,
				new volScalarField
				(		       
					IOobject
					(
					   "fluxYMo"+Foam::name(fluxYI)+Foam::name(l)+Foam::name(r),
					    runTime.timeName(),
					    mesh,
					    IOobject::READ_IF_PRESENT,
					    IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
        }
}

//Create average angular neutron fluxes
Info<< "Creating average angular flux\n" << endl;
PtrList<PtrList<PtrList<volScalarField> > > psi(energyGroups);
forAll(psi, energyI)
{
	psi.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(psi[energyI], i)
	{
		psi[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(psi[energyI][i], j)
		{
			psi[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"psi"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-2,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

PtrList<PtrList<PtrList<volScalarField> > > psiX(energyGroups);
forAll(psiX, energyI)
{
	psiX.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(psiX[energyI], i)
	{
		psiX[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(psiX[energyI][i], j)
		{
			psiX[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"psiX"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

PtrList<PtrList<PtrList<volScalarField> > > psiY(energyGroups);
forAll(psiY, energyI)
{
	psiY.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(psiY[energyI], i)
	{
		psiY[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(psiY[energyI][i], j)
		{
			psiY[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"psiY"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

PtrList<PtrList<PtrList<volScalarField> > > psiH(energyGroups);
forAll(psiH, energyI)
{
	psiH.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(psiH[energyI], i)
	{
		psiH[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(psiH[energyI][i], j)
		{
			psiH[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"psiH"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}
						
//Create total neutron source
Info<< "Total isotropic neutron source\n" << endl;
PtrList <volScalarField> Q(energyGroups);
forAll(Q, i)
{
	Q.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Q"+Foam::name(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

//Create angular neutron source for anisotropic scattering
Info<< "Angular neutron source\n" << endl;
PtrList<PtrList<PtrList <volScalarField> > > qm(energyGroups);
forAll(qm, energyI)
{
	qm.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(qm[energyI], i)
	{
		qm[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(qm[energyI][i], j)
		{
			qm[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"qm"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}
					

//Create neutron source x-moment
Info<< "Total neutron source x-moment\n" << endl;
PtrList <volScalarField> QX(energyGroups);
forAll(QX, i)
{
	QX.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"QX"+Foam::name(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-4,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

PtrList <volScalarField> qX(energyGroups);
forAll(qX, i)
{
	qX.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"qX"+Foam::name(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-4,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

Info<< "Total neutron source y-moment\n" << endl;
PtrList <volScalarField> QY(energyGroups);
forAll(QY, i)
{
	QY.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"QY"+Foam::name(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-4,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

PtrList <volScalarField> qY(energyGroups);
forAll(qY, i)
{
	qY.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"qY"+Foam::name(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-4,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

//Create angular neutron source for anisotropic scattering
Info<< "Angular neutron source x-moment\n" << endl;
PtrList<PtrList<PtrList <volScalarField> > > QmX(energyGroups);
forAll(QmX, energyI)
{
	QmX.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(QmX[energyI], i)
	{
		QmX[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(QmX[energyI][i], j)
		{
			QmX[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"QmX"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-4,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

PtrList<PtrList<PtrList <volScalarField> > >qmX(energyGroups);
forAll(qmX, energyI)
{
	qmX.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(qmX[energyI], i)
	{
		qmX[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(qmX[energyI][i], j)
		{
			qmX[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"qmX"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-4,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

Info<< "Angular neutron source y-moment\n" << endl;
PtrList<PtrList<PtrList <volScalarField> > > QmY(energyGroups);
forAll(QmY, energyI)
{
	QmY.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(QmY[energyI], i)
	{
		QmY[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(QmY[energyI][i], j)
		{
			QmY[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"QmY"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-4,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

PtrList<PtrList<PtrList <volScalarField> > > qmY(energyGroups);
forAll(qmY, energyI)
{
	qmY.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(qmY[energyI], i)
	{
		qmY[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(qmY[energyI][i], j)
		{
			qmY[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"qmY"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-4,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}



Info<< "fissionSource\n" << endl;
PtrList<volScalarField> fissionSource(energyGroups);
forAll(fissionSource, i)
{
	fissionSource.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"fissionSource"+Foam::name(i),
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

Info<< "totalFissions\n" << endl;
volScalarField totalFissions
(
	IOobject
	(
		"totalFissions",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);

scalar volFissions=0.0;
scalar prevFissions=volFissions;

Info<< "scatteringSource\n" << endl;
PtrList<volScalarField> scatteringSource(energyGroups);
forAll(scatteringSource, i)
{
	scatteringSource.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"scatteringSource"+i,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

PtrList <volScalarField> oldQ(energyGroups);
PtrList <volScalarField> prevFlux(energyGroups);
PtrList <volScalarField> oldFlux(energyGroups);

//Define fields for C factors for linear source calculations
Info<< "C factors\n" << endl;
PtrList<volScalarField> Cxx(energyGroups);
forAll(Cxx, i)
{
	Cxx.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Cxx"+i,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

PtrList<volScalarField> Cyy(energyGroups);
forAll(Cyy, i)
{
	Cyy.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Cyy"+i,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

PtrList<volScalarField> Cxy(energyGroups);
forAll(Cxy, i)
{
	Cxy.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"Cxy"+i,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

//  CHECK UNITS!! //
//Define fields for C factors for anisotropic linear source calculations
Info<< "Cm factors\n" << endl;
PtrList<PtrList<PtrList <volScalarField> > > Cmy(energyGroups);
forAll(Cmy, energyI)
{
	Cmy.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(Cmy[energyI], i)
	{
		Cmy[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(Cmy[energyI][i], j)
		{
			Cmy[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"Cmy"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-2,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

PtrList<PtrList<PtrList <volScalarField> > > Cmx(energyGroups);
forAll(Cmx, energyI)
{
	Cmx.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(Cmx[energyI], i)
	{
		Cmx[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(Cmx[energyI][i], j)
		{
			Cmx[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"Cmx"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-2,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

PtrList<PtrList<PtrList <volScalarField> > > Cm(energyGroups);
forAll(Cm, energyI)
{
	Cm.set(energyI, new PtrList<PtrList<volScalarField> >(naz));
	forAll(Cm[energyI], i)
	{
		Cm[energyI].set(i, new PtrList<volScalarField>(npo));
		forAll(Cm[energyI][i], j)
		{
			Cm[energyI][i].set
			(
				j,
				new volScalarField
				(
					IOobject
					(
						"Cm"+Foam::name(energyI)+Foam::name(i)+Foam::name(j),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-2,-1,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

//Define factors for calculating QX from qX and qY etc
Info<< "M factors\n" << endl;
volScalarField Mxx
(
	IOobject
	(
		"Mxx",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);

volScalarField Myy
(
	IOobject
	(
		"Myy",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);

volScalarField Mxy
(
	IOobject
	(
		"Mxy",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);

Info<< "Nuclear data fields\n" << endl;
Info<< "nuSigmaEff\n" << endl;
PtrList<volScalarField> nuSigmaEff(energyGroups);
forAll(nuSigmaEff, i)
{
	nuSigmaEff.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"nuSigmaEff"+i,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

Info<< "sigmaT\n" << endl;
PtrList<volScalarField> sigmaT(energyGroups);
forAll(sigmaT, i)
{
	sigmaT.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"sigmaT"+i,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

Info<< "sigmaA\n" << endl;
PtrList<volScalarField> sigmaA(energyGroups);
forAll(sigmaA, i)
{
	sigmaA.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"sigmaA"+i,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

Info<<"sigmaS\n" << endl;
PtrList<PtrList<PtrList<volScalarField> > > sigmaS(energyGroups);
forAll(sigmaS, energyI)
{
	sigmaS.set(energyI, new PtrList<PtrList<volScalarField> >(energyGroups));
	forAll(sigmaS[energyI], energyJ)
	{
		sigmaS[energyI].set(energyJ, new PtrList<volScalarField>(anisotropy+1));
		forAll(sigmaS[energyI][energyJ], l)
		{
			sigmaS[energyI][energyJ].set(
				l,
				new volScalarField
				(
					IOobject
					(
						"sigmaS"+Foam::name(energyI)+Foam::name(energyJ)+Foam::name(l),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("", dimensionSet(0,-1,0,0,0,0,0), 0.0),
					zeroGradientFvPatchScalarField::typeName
				)
			);
		}
	}
}

Info<< "chi\n" << endl;
PtrList<volScalarField> chi(energyGroups);
forAll(chi, i)
{
	chi.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"chi"+i,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}



//Read geometry and create constants
Info<< "Reading width, height and depth of geometry"<<endl;
boundBox cuboidDomain(mesh.points());
point extreme=cuboidDomain.max();
scalar w=extreme.x();
scalar h=extreme.y();
scalar depth=extreme.z();
scalar d_over_2=depth/2.0;

scalarField vol=mesh.V();
scalarField area=vol/depth;

//For calculating cell centres and moments in linear source
volVectorField centres=mesh.C();
