    ///MRF.correctBoundaryVelocity(U);
  
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      - fvm::Sp(contErr, U)
      + MRF.DDt(rho, U)//MRF是多参考坐标系，在动网格中会用到（像旋转机械等），一般情况下可忽略
      + turbulence.divDevRhoReff(U)
      + fvm::Sp(rho*coupling.damping(), U)    /////
      + fvm::Su(rho*coupling.source(), U)     /////     
      
     ==
        fvOptions(rho, U)//源项或约束
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
