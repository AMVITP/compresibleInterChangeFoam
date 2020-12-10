
//-------------------------------------------------------------------------------
volScalarField TTsat
(
    IOobject
    (
        "TTsat",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    1/(1/Tsat -RRv/hEvap*log(p/Psat))
);
//--------------------------------------------------------------------------------

