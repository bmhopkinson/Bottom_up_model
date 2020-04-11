function [ P ] = PvE( E, Params )
%paturating exponential function to calculate photosynthesis as a function of irradiance
    Pmax = Params(1);
    Ek   = Params(2);
    R    = Params(3);

    P = Pmax.*(1-exp(-E./Ek)) - R;

end

