function [ zs,Ints,dAs,dTs,Rtotal,Ttotal ] = multilayerAbsorption(layers, incidence, wavelength, steps )
%multilayerAbsorption Calculates the intensity, absorption and temperature increase profiles 
%in each layer of a multilayers structure for p-polarized light.
%  Arguments:
%        layers: a 2D array, each rows corresponds to one layer, starting from
%            the top layer, usually vacuum or air, down to the bottom layer,
%            usually the substrate. For each layer, i.e. row, the complex index
%            of refraction n+jk, the density [g.cm-3], the heat capacity
%            [J.kg-1.K-1] and the thickness [m] are given in that order.
%        incidence: the angle of incidence, in radians.
%        wavelength: the vacuum wavelength [m].
%        steps: the number of data points for the calculation of the intensity,
%            absorption, temperature increase profiles for each layer. If steps
%            is None, then only Rtotal and Ttotal are calculated
    
%    Returns:
%        zs: a 2D array containing the depth z [m] within each layer where the
%            other quantities are calculated.
%        Ints: a 2D array containing the intensity profiles within each layer
%            for an incoming intensity is 1 J.m-2 or 0.1 mJ.cm-2
%        dAs: a 2D array containing the differential absorption [m-1] within
%            each layer.
%        dTs: a 2D array containing the temperature increase [K] within each
%            layer.
%        Rtotal: total amount of reflection from the multilayer.
%        Ttotal: total transmission in the last layer of the multilayer.

 nblayers = length(layers(:,1));
        N = zeros(nblayers,1);
thickness = zeros(nblayers,1);
  density = zeros(nblayers,1);
heatcapacity = zeros(nblayers,1);

for n = 1:nblayers
    N(n)                  = layers(n,1) + 1i*layers(n,2);               %refractive index n+ik
    density(n)            = 1e3*layers(n,3);                            %density in g/cm3 or kg/m3
    heatcapacity(n)       = layers(n,4);                                %heat capacity in J/kg.K
    thickness(n)          = layers(n,5);                                %thickness [m]
end%for
    
%Snell laws
theta    =  zeros(nblayers, 1);
theta(1) =  incidence;

for n = 2:nblayers
    theta(n) = asin(N(1)/N(n)*sin(theta(1)));
end%for

%fresnel coefficient for P polarized light
rfresnel = zeros(nblayers -1,1);
tfresnel = zeros(nblayers -1,1);

for n = 1:(nblayers-1)
    rfresnel(n) = (N(n+1)*cos(theta(n)) - N(n)*cos(theta(n+1)))/(N(n+1)*cos(theta(n)) + N(n)*cos(theta(n+1)));
    tfresnel(n) = 2.0*N(n)*cos(theta(n))/(N(n+1)*cos(theta(n)) + N(n)*cos(theta(n+1)));
end%for

%interface change matrix
Jnm = zeros(2,2,nblayers-1);
for n = 1:nblayers-1
    Jnm(1,1,n) = 1.0/tfresnel(n);
    Jnm(1,2,n) = rfresnel(n)/tfresnel(n);
    Jnm(2,1,n) = rfresnel(n)/tfresnel(n);
    Jnm(2,2,n) = 1/tfresnel(n);
end;%for

%calculating z-component of the wave vector
Kz = 2.0.*pi./(wavelength).*N.*cos(theta);

% #phase changes
beta      = Kz.*thickness;
Ln        = zeros(2,2,nblayers-1);
Ln(:,:,1) = [1,0 ; 0,1];

for n = 2:nblayers-1
     Ln(1,1,n) = exp(-1i*beta(n));
     Ln(1,2,n) = 0;
     Ln(2,1,n) = 0;
     Ln(2,2,n) = exp(1i*beta(n));
end;%for

% calculating propagation matrix
S = Jnm(:,:,nblayers-1);
for n =  (nblayers-2):-1:1
    S = Jnm(:,:,n)*(Ln(:,:,n+1)*S);
end%for

%Total transmission and reflection of the multilayer
Rtotal = abs(S(2,1)/S(1,1))^2;
Ttotal = real(conj(N(nblayers))*cos(theta(nblayers))/(N(1)*cos(theta(1))))*abs(1/S(1,1))^2;
%%

%calculating D matrix for intermediate field
Dn = zeros(2,2,nblayers);
Dn(1,1,nblayers) = 1.0/S(1,1);
Dn(1,2,nblayers) = 0.0;
Dn(2,1,nblayers) = 0.0;
Dn(1,1,nblayers) = 1.0/S(1,1);

for n =  (nblayers-1):-1:1
    Dn(:,:,n) = Ln(:,:,n)*(Jnm(:,:,n)*Dn(:,:,n+1));
end;%for  

%For each layer except the first, compute Intensity,Absorption, Temperature
zs   = zeros(nblayers, steps);
Ints = zeros(nblayers, steps);
dAs  = zeros(nblayers, steps);
dTs  = zeros(nblayers, steps);

for n = 1:(nblayers)
    zs(n,:) = logspace(-10.0, log10(thickness(n)), steps);
    Ep      = Dn(1,1,n)*exp(1i*Kz(n)*zs(n,:));
    Em      = Dn(2,1,n)*exp(-1i*Kz(n)*zs(n,:));
    Etx     = 0.0*Ep + 0.0*Em;
    Ety     =  cos(theta(n))*Ep - cos(theta(n))*Em;
    Etz     = -sin(theta(n))*Ep - sin(theta(n))*Em;
    Ints(n,:) = real(Etx .* conj(Etx) + Ety.*conj(Ety) + Etz.*conj(Etz));
    dAs(n,:)  = real(N(n)*cos(theta(n))/(N(1)*cos(theta(1)))*2*imag(Kz(n))*Ints(n,:));
    dTs(n,:)  = dAs(n,:)/(density(n)*heatcapacity(n));
end;%for


    
end

