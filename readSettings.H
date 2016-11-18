Info<<"Reading angles, boundary albedos and settings - MOCsettings\n"<<endl;
IOdictionary MOCsettings
(
    IOobject
    (
        "MOCsettings",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//Read line spacing
scalar ts (readScalar(MOCsettings.lookup("spacing")));       

 //Read azimuthal angles
int naz (readInt(MOCsettings.lookup("naz")));                 

//Read polar angles
const int npo(readInt(MOCsettings.lookup("npo")));             

//Read left boundary albedo
const scalar alphaL(readScalar(MOCsettings.lookup("alphaL"))); 

//Read right boundary albedo
const scalar alphaR(readScalar(MOCsettings.lookup("alphaR"))); 

//Read top boundary albedo
const scalar alphaT(readScalar(MOCsettings.lookup("alphaT"))); 

//Read bottom boundary albedo
const scalar alphaB(readScalar(MOCsettings.lookup("alphaB"))); 

//Read error tolerance
const scalar tol(readScalar(MOCsettings.lookup("tol"))); 

//Read whether to use linear source
const bool linearSource=MOCsettings.lookupOrDefault("linearSource",0);

//Read whether to integrate zone fluxes
const bool zoneIntegral=MOCsettings.lookupOrDefault("zoneIntegral",0);


//Read order of anisotropy
/*
int anisotropy(readInt(MOCsettings.lookupOrDefault("anisotropy",0)))
if((anisotropy<0) || (anisotropy>5))
{
	Info<<"Order of anisotropy must be between 0 and 5 - reset to 0"<<endl;
	anisotropy=0;
}
*/
