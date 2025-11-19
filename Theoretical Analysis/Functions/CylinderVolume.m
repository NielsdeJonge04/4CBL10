function V = CylinderVolume (Ca , Cyl )
    % Geometry
    B = Cyl.Bore ;
    S = Cyl.Stroke ;
    CR = Cyl.CompressionRatio ;
    Lcr = Cyl.ConRod ;
    Ap = pi*(B^2)/4; % piston area
    Vd = Ap * S ; % displacement volume
    Vc = Vd / ( CR - 1) ; % clearance volume

    % Shift angle so that TDC = 0
    theta = deg2rad(Ca - Cyl.TDCangle ) ;

    % Slider crank kinematics
    r = S /2;
    under = Lcr^2 - (r*sin(theta)).^2;
    under ( under < 0) = 0;
    x = r *(1-cos(theta))+(Lcr - sqrt(under));

    % Total cylinder volume
    V = Vc + Ap * x ;
end