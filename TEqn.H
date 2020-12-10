//opn2

//volScalarField k = alpha1*k1+(scalar(1)-alpha1)*k2;

// turbulent heat transfer is to be added some time

//volScalarField rhoCp = alpha1*rho1*cp1+(scalar(1)-alpha1)*rho2*cp2;

volScalarField Cp1= mixture.thermo1().Cp();
volScalarField Cp2= mixture.thermo2().Cp();


volScalarField rhoCp = alpha1*rho1*Cp1+alpha2*rho2*Cp2;

//volScalarField CpAlph = alpha1*Cp1+alpha2*Cp2;

//volScalarField CpAlph1 = alpha1/Cp1+alpha2/Cp2;
//volScalarField CpAlph = 1/CpAlph1;


volScalarField CpA = rhoCp/rho;


rhoCp.oldTime() = alpha1.oldTime()*rho1*Cp1+alpha2.oldTime()*rho2*Cp2;


//rhoCp.oldTime() = alpha1.oldTime()*rho1*cp1+(scalar(1)-alpha1.oldTime())*rho2*cp2;

//ins ----------------------------------------------------------------------------
//volScalarField CpAlph = alpha1*cp1+(scalar(1)-alpha1)*cp2;
//CpAlph.oldTime() = alpha1.oldTime()*cp1+(scalar(1)-alpha1.oldTime())*cp2;
//volScalarField phiCp = alpha1*phi1*cp1+(scalar(1)-alpha1)*phi2*cp2;
//phiCp.oldTime() = alpha1.oldTime()*phi1*cp1+(scalar(1)-alpha1.oldTime())*phi2*cp2;



//-------------------------------------------------------------------------------

//surfaceScalarField alphaPhi = (rhoPhi-phi*rho2)/(rho1-rho2);
//surfaceScalarField rhoCpPhi = alphaPhi1*rho1*Cp1+alphaPhi2*rho2*Cp2;

// surfaceScalarField rhoCpPhi = alphaPhi*(rho1*cp1-rho2*cp2)+phi*rho2*cp2;

//ins ---------------------------------------------------------------------------
//surfaceScalarField rhoCpPhi( "rhoCpPhi",   fvc::interpolate(CpAlph)*rhoPhi);

surfaceScalarField rhoCpPhi( "rhoCpPhi",   fvc::interpolate(rhoCp)*phi);


//surfaceScalarField ktutbCpAlph(  "ktutbCpAlph", fvc::interpolate(CpAlph)*fvc::interpolate(mixture.alphaEff(turbulence->mut())));

surfaceScalarField ktutbCpAlph(  "ktutbCpAlph", fvc::interpolate(CpA)*fvc::interpolate(mixture.alphaEff(turbulence->mut())));



//volScalarField::Internal contErrCpAlph("contErrCpAlph", CpAlph*contErr);
volScalarField::Internal contErrCpAlph("contErrCpAlph", CpA*contErr);



//surfaceScalarField absphiUCpAlph("absphiUCpAlph", fvc::interpolate(CpAlph)*fvc::absolute(phi, U));

surfaceScalarField absphiUCpAlph("absphiUCpAlph", fvc::interpolate(CpA)*fvc::absolute(phi, U));

//-------------------------------------------------------------------------------
//opn2
// source terms due to mass sources (due to non divergence-free velocity field)
//hCSource = fvc::ddt(rhoCp) + fvc::div(rhoCpPhi);

// energy equation
fvScalarMatrix TEqn
(
        fvm::ddt(rhoCp,T) + fvm::div(rhoCpPhi,T) 
//opn2 -------------------------------------------------------
//      - fvm::laplacian(k,T) 

// opn5 -------------------------------------------------------
//     - fvm::Sp(contErr, T)
//   - fvm::laplacian(mixture.alphaEff(turbulence->mut()), T)
//ins
      - fvm::Sp(contErrCpAlph, T)
     - fvm::laplacian(ktutbCpAlph, T)

// -----------------------------------------------------------
      - fvm::Sp(hESource,T) 
//      - fvm::Sp(hCSource,T) 
//opn5 -------------------------------------------------------

     +  (
            //fvc::div(fvc::absolute(phi, U), p)
              fvc::div(absphiUCpAlph, p)
            + fvc::ddt(rhoCp, K) + fvc::div(rhoCpPhi, K)
        )
       *(
           alpha1/mixture.thermo1().Cv()
         + alpha2/mixture.thermo2().Cv()
        )

//-----------------------------------------------------------
      == -hESource*TTsat
);

TEqn.solve();






// opn5
/*
{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T)
      - fvm::Sp(contErr, T)
      - fvm::laplacian(mixture.alphaEff(turbulence->mut()), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
        )
       *(
           alpha1/mixture.thermo1().Cv()
         + alpha2/mixture.thermo2().Cv()
        )
     ==
        fvOptions(rho, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    mixture.correctThermo();
    mixture.correct();
}
*/
