forAll(QX, energyI)
{
	QX[energyI] *= 0;
	QY[energyI] *= 0;

	forAll(fluxX,energyJ)
	{
		//Calculate source moments		
		QX[energyI] += (sigmaS[energyI][energyJ][0] + chi[energyI] * nuSigmaEff[energyJ]/keff) * fluxX[energyJ];
		QY[energyI] += (sigmaS[energyI][energyJ][0] + chi[energyI] * nuSigmaEff[energyJ]/keff) * fluxY[energyJ];
	}

	//Calculate x and y expansion coefficients by inverting the governing linear system:
	//M*qxy = Qxy
	qX[energyI].internalField() = (Myy.internalField()*QX[energyI].internalField()-Mxy.internalField()*QY[energyI].internalField())
				/(Myy.internalField()*Mxx.internalField() -Mxy.internalField()*Mxy.internalField());
	qY[energyI].internalField() = (Mxx.internalField()*QY[energyI].internalField()-Mxy.internalField()*QX[energyI].internalField())
				/(Myy.internalField()*Mxx.internalField() -Mxy.internalField()*Mxy.internalField());
}

