
//	volScalarField Pr = Prc*rho;//Prc*alphac + Prd*(1.0 - alphac);
//	volScalarField Prt = Prtc*rho;//Prtc*alphac + Prtd*(1.0 - alphac);

      kappat = rho* turbulence->nut()*Cp/Prt;
      kappat.correctBoundaryConditions();

    volScalarField kappaEff("kappaEff", K + kappat);

    fvScalarMatrix TEqn
    (
           fvm::ddt(rho*Cp,T) 
         //- fvm::Sp(fvc::ddt(rho),T)
         + fvm::div(fvc::interpolate(Cp)*rhoPhi,T) 
         //- fvm::Sp(fvc::interpolate(rho)*phi,T)
         - fvm::laplacian(kappaEff,T)
    );

    TEqn.relax();
    TEqn.solve();
     Cp =  0.833 +(4.141 * T * dimSet1);
     //K  =  (0.163 - 1.244 * 1.0e-04 * T * dimSet1)*dimSet2;
    //K=(-0.5744272+(5.948806 * T * dimSet1*1.0e-03)-(1.66063 * pow(T,2) * dimSet1* dimSet1*1.0e-05)+(1.507483 * pow(T,3) * dimSet1* dimSet1*dimSet1*1.0e-08))*dimSet2;
    //rho  =  (820.504 + 1.23 * 1.0e-01 * T * dimSet1);
    //mu=(0.956994130219100180000-(9.34196070131963940000 * T * dimSet1*1.0e-03)+(3.43772166366523290000 * pow(T,2) * dimSet1* dimSet1*1.0e-05)-(5.64114506183098770000 * pow(T,3) * dimSet1* dimSet1*dimSet1*1.0e-08)+(3.47903390933215070000 * pow(T,4) *dimSet1* dimSet1* dimSet1*dimSet1*1.0e-011))*dimSet2;
//mu=(1.026608-(0.01011998 * T * dimSet1)+(3.760667 * pow(T,2) * dimSet1* dimSet1*1.0e-05)-(6.234743 * pow(T,3) * dimSet1* dimSet1*dimSet1*1.0e-08)+(3.886479 * pow(T,4) *dimSet1* dimSet1* dimSet1*dimSet1*1.0e-011))*dimSet2;
    //Info<< "CCCCCCCCCCCCpppp"<<Cp<<endl;
    Info<< "Temperature Field results = "
           << T.weightedAverage(mesh.V()).value()
           << "  Min(T) = " << min(T).value()
           << "  Max(T) = " << max(T).value()
           << endl;

