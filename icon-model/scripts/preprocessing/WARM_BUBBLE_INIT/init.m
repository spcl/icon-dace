% ICON
%
% ------------------------------------------
% Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
% Contact information: icon-model.org
% See AUTHORS.TXT for a list of authors
% See LICENSES/ for license information
% SPDX-License-Identifier: BSD-3-Clause
% ------------------------------------------

%initialize warm bubble case
%Ref: G. H. Bryan and J. M. Fritsch, "A benchmark simulation for moist
%nonhydrostatic numerical models", MWR, 2002

clear; clc;

%Constants
cpd = 1004; cpl = 4186; cpv = 1885;
g = 9.81; Rd = 287; Rv = 461; e = Rd/Rv;
rd_o_cpd = Rd/cpd;
pref = 1e5; 

ps = pref; pt = 1000; %surface and top pressure in Pa

the = 320;
rt  = 0.02;

p = ps:-100:pt; %in pa
n = length(p);
temp0 = 300;

th = zeros(n,1);
rv = zeros(n,1);
pi = zeros(n,1);

for i = 1:n
    
    fun = @(temp) findzero(the,p(i),temp,cpl,cpv,e,rt,cpd,Rd,pref);
    
    temp = fzero(fun,temp0);       
    
    pi(i)  = (p(i)/pref)^rd_o_cpd;
    th(i)  = temp/pi(i);    
    esat = 610.78 * exp( (17.2694 * (temp-273.15))/(temp-35.86) );
    rv(i) = 0.622 * esat / (p(i)-0.378*esat);
   
    temp0 = temp;
end;

%Get cloud liquid water
rc = rt-rv;

%Get height
thr = th .* (1 + rv/e)/(1+rt);
z = zeros(n-1,1);
for i = 2:n
    z(i) = z(i-1) - (cpd/g)*(thr(i)+thr(i-1))*0.5*(pi(i)-pi(i-1));
end;

data = [z th rv rc zeros(n,1)];
str = strcat('sound_WarmBubble_the_',num2str(the),'_rt_',num2str(rt));
save(str,'data','-ascii');

