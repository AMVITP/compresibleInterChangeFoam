#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"
//-----------------------------------------------------------------------------------
// opn2 

   //- Create fields needed for phase change model
    Info<< "Reading field psi0\n" << endl;
    volScalarField psi0
    (
        IOobject
        (
            "psi0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        psi0
    );

        volScalarField rhoSource
        (
                IOobject
                (
                        "rhoSource",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("rhoSource", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
        );

       volScalarField hESource
        (
                IOobject
                (
                        "hESource",
                        runTime.timeName(),
                        mesh
                ),
                mesh,
                dimensionedScalar("hESource", dimensionSet(1,-1,-3,-1,0,0,0), 0.0)
        );

        volScalarField hCSource
        (
                IOobject
                (
                        "hCSource",
                        runTime.timeName(),
                        mesh
                ),
                mesh,
                dimensionedScalar("hCSource", dimensionSet(1,-1,-3,-1,0,0,0), 0.0)
        );

        volScalarField alphaSourceSu
        (
                IOobject
                (
                        "alphaSourceSu",
                        runTime.timeName(),
                        mesh
                ),
                mesh,
                dimensionedScalar("alphaSourceSu", dimensionSet(0,0,-1,0,0,0,0), 0.0)
        );

        volScalarField alphaSourceSp
        (
                IOobject
                (
                        "alphaSourceSp",
                        runTime.timeName(),
                        mesh
                ),
                mesh,
                dimensionedScalar("alphaSourceSp", dimensionSet(0,0,-1,0,0,0,0), 0.0)
        );

        dimensionedScalar mEvap ("mEvap", dimensionSet(1,0,-1,0,0,0,0), 0.0);

//------------------------------------------------------------------------------------

















Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Reading thermophysical properties\n" << endl;

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();
//----------------------------------------------------------------------------------
//

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

        dictionary phase1 = transportProperties.subDict("phase1");
        dictionary phase2 = transportProperties.subDict("phase2");
        dimensionedScalar k1 = phase1.lookup("k");
        dimensionedScalar cp1 = phase1.lookup("cp");
        dimensionedScalar k2 = phase2.lookup("k");
        dimensionedScalar cp2 = phase2.lookup("cp");
        dimensionedScalar Rph = transportProperties.lookup("Rph");
        dimensionedScalar Tsat = transportProperties.lookup("Tsat");
        dimensionedScalar hEvap = transportProperties.lookup("hEvap");
        dimensionedScalar DPsi = transportProperties.lookup("DPsi");
        dimensionedScalar DAlpha = transportProperties.lookup("DAlpha");
        dimensionedScalar Psat = transportProperties.lookup("Psat");
        dimensionedScalar RRv = transportProperties.lookup("RRv");





//--------------------------------------------------------------------------------
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
//--------------------------------------------------------------------------------
  rho.oldTime();
//-------------------------------------------------------------------------------

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"







// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt(alpha1*fvc::div(phi));











// Construct compressible turbulence model
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New(rho, U, rhoPhi, mixture)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"
