volScalarField hEvapP=hEvap*p;

//volScalerField RRvTTsat=RRv*TTsat;

//volScalerField RRvTTsat2=RRvTTsat*TTsat;

solve (   fvm::ddt(hEvapP,TTsat) == RRv*TTsat*TTsat*fvc::ddt(p));
