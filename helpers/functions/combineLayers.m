function [z, Int, dA, dT] = combineLayers(zs, Ints, dAs, dTs)
%Combine all the 2D array into a single 1D array.
%
% Arguments:
% zs:   a 2D array containing the depth within each layer.
% Ints: a 2D array containing the intensity profiles within each layer.
% dAs:  a 2D array containing the differential absorption within each layer.
% dTs:  a 2D array containing the temperature increase within each layer.
%     
% Returns:
% z: a 1D array containing the depth within the whole multilayer structure.
% Int: a 1D array containing the intensity profiles within the whole multilayer structure.
% dA: a 1D array containing the differential absorption within the whole multilayer structure.
% dT: a 1D array containing the temperature increase within the whole multilayer structure.

nblayers =    length(zs(:,1));
steps    =    length(zs(1,:));
nzs      =               zs;
for k = 2:1:nblayers
    %disp(num2str(k));
    nzs(k,:) = nzs(k,:) + nzs(k-1,steps);
end

z   = reshape(nzs',  nblayers*steps,1);
Int = reshape(Ints', nblayers*steps,1);
dA  = reshape(dAs',  nblayers*steps,1);
dT  = reshape(dTs',  nblayers*steps,1);

end

