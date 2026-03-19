
    Info<< "Reading transportProperties\n" << endl;
    IOdictionary transportPropertiesT
    (
        IOobject
        (
            "transportPropertiesT",
            runTime.constant(),
            newMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
       Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            newMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        newMesh
    );
    Info<< "Reading field gammaT\n" << endl;
       volScalarField gammaT
    (
        IOobject
        (
            "gammaT",
            runTime.timeName(),
            newMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        newMesh,
        dimensionedScalar("one", dimless, 1)
    );
    dimensionedScalar DTS //add
    (
        transportPropertiesT.lookup("DTS")
    );   
    dimensionedScalar DTL //add
    (
        transportPropertiesT.lookup("DTL")
    ); 
    volScalarField DDD
    (
        IOobject
        (
            "DDD",
            runTime.timeName(),
            newMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        newMesh,
        dimensionedScalar("one", dimArea/dimTime, 1)
    ); 
        surfaceScalarField phiT
        (
            IOobject
            (
                "phiT",
                runTime.timeName(),
                newMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            newMesh,
            dimensionedScalar("zero", dimArea*dimVelocity, 0),
            calculatedFvsPatchScalarField::typeName
        );
    Info<< "Reading field gammaT\n" << endl;
    volScalarField Trev
    (
        IOobject
        (
            "Trev",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimTemperature, 1)
    );
    volScalarField TIbFace
    (
        IOobject
        (
            "TIbFace",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTemperature, 0)
    );
    dimensionedScalar TRef //add
    (
        transportPropertiesT.lookup("TRef")
    );  


    volScalarField gammas
    (
        IOobject
        (
            "gammas",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimTemperature, 1)
    );
//------------------------------------------------

