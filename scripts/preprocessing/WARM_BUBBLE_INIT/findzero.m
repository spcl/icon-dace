% ICON
%
% ------------------------------------------
% Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
% Contact information: icon-model.org
% See AUTHORS.TXT for a list of authors
% See LICENSES/ for license information
% SPDX-License-Identifier: BSD-3-Clause
% ------------------------------------------

%saturation mixing ratio as a function of T and P

function y = findzero(the,pres,temp,cpl,cpv,e,rt,cpd,Rd,pref)

    y = the - temp*reverse_pid(temp,pres)*exp(latent_heat(temp)*sat_humidity(temp,pres)/((cpd+cpl*rt)*temp));

    function Lv = latent_heat(temp)
        Lv = 2.5e6 - (cpl - cpv)*(temp-273.15);
    end
    
    function rpid = reverse_pid(temp,pres)
        pd    = pres*e/(sat_humidity(temp,pres)+e);
        rpid  = (pd/pref)^(-Rd/(cpd+cpl*rt)); %reverse pid
    end
    
    function rvs = sat_humidity(temp,pres)      
        esat = 610.78 * exp( (17.2694 * (temp-273.15))/(temp-35.86) );
        rvs = 0.622 * esat / (pres-0.378*esat);
    end
end