//------------------------------------------------------------------------------------
// opn 2

/*opn2
    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    volScalarField muEffC
    (
        "muEffC",

        twoPhaseProperties.mu()
      + rho*turbulence->nut()
    );

*/
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
//opn5
     - fvm::Sp(contErr, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
//

/*opn2
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    //  additional term because velocity field is not divergence-free

        - muEffC*fvc::grad(fvc::div(U)) 
*/

    );

    UEqn.relax();


// opn2
//    if (momentumPredictor)

// opn5
    if (pimple.momentumPredictor())
// 

    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
//                  fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)

// opn2             fvc::interpolate(sigmaK)*fvc::snGrad(alpha1)
                    mixture.surfaceTensionForce()        
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }




//------------------------------------------------------------------------------------





//------------------------------------------------------------------------------------
//  opn5 
/*
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(contErr, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);

        K = 0.5*magSqr(U);
    }
*/
//----------------------------------------------------------------------------------
