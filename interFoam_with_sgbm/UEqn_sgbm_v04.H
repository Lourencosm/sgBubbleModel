    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        if (bubbleModel_->isActive())              // LSM add: air_water entrainment
        {    
            // buoyancy force v0
//            const volScalarField& Cb_ = mesh.lookupObject<volScalarField>("Cb");
//            volVectorField buoyancyForce = rho*Cb_*(-g);

            // buoyancy force v1
//            const volScalarField& Cb_bounded_ = mesh.lookupObject<volScalarField>("Cb_bounded"); 
//            volVectorField buoyancyForce = rho*Cb_bounded_*(-g);

//            // buoyancy force v2
//            const volScalarField& alphaWater_ = mesh.lookupObject<volScalarField>("alpha.water"); 
//            const volScalarField& Cb_bounded_ = mesh.lookupObject<volScalarField>("Cb_bounded"); 
//            volVectorField buoyancyForce = rho*alphaWater_*Cb_bounded_*(-g);

            // buoyancy force v3
            const volScalarField& alphaWater_ = mesh.lookupObject<volScalarField>("alpha.water"); 
            const volScalarField& Cb_ = mesh.lookupObject<volScalarField>("Cb"); 
            volVectorField buoyancyForce = rho*alphaWater_*Cb_*(-g);

            Info<< "    buoyancyForce - min, max, avg       : " << gMin(buoyancyForce) << " " << gMax(buoyancyForce) << " " << gAverage(buoyancyForce) << endl;
            //    Info<< "    buoyancyForce - dimensions          : " << buoyancyForce.dimensions() << endl;

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
                + buoyancyForce
            );
        }
        else
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
        }



        fvOptions.correct(U);
    }
