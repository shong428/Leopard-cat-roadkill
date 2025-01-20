function [Bmus,qerror_geo,qerror_dta] = som_bmus_geo2(sMap, sData, kradius)
% function [Bmus,qerror_geo,qerror_dta] = som_bmus_geo2(sMap, sData, kradius)
%
% Find the BMU for each data pattern, together with the geographical error,
% and the error in the NON-Geographic variables
%
% INPUTS
%   sMAP    - SOM map
%   sData   - SOM data
%   kradius - Geographic tolerance (squared)
%   
% OUTPUTS
%   Bmus        - Index of the BMU for each pattern
%   qerror_geo  - Geographic quantization error (BMU-Pattern)
%   qerror_dta  - Data (NON-Geographic) quantization error (BMU-Pattern)

% V1.0 by V.Lobo & F.Bacao, 31/10/2004


% SOM_BMUS Find the best-matching units from the map for the given vectors.
% check arguments and initialize

% ---- GEO initializations
geo=0;Ud = som_unit_dists(sMap.topol).^2;
% kradius2=kradius^2;
% ---- End GEO initializations;
error(nargchk(1, 3, nargin));  % check no. of input args is correct

if kradius == -1
    kradius = inf;
end;

% sMap
if isstruct(sMap), 
  switch sMap.type, 
   case 'som_map', M = sMap.codebook; 
   case 'som_data', M = sMap.data;
   otherwise, error('Invalid 1st argument.');
  end
else 
  M = sMap; 
end
[munits dim] = size(M);
if any(any(isnan(M))), 
  error ('Map codebook must not have missing components.');
end

% data
if isstruct(sData), 
  switch sData.type, 
   case 'som_map', D = sData.codebook;
   case 'som_data', D = sData.data;
   otherwise, error('Invalid 2nd argument.');
  end
else 
  D = sData;
end
[dlen ddim] = size(D);
if dim ~= ddim, 
  error('Data and map dimensions do not match.')
end

qerror_geo=zeros(1,dlen);
qerror_dta=qerror_geo;
Bmus=qerror_geo;

for i=1:dlen        % search each data point independently
    Dst=dist( D(i,1:2),M(:,1:2)');
    [qerror_geo(i) bmu_geo]=min(Dst);

    non_searchable_units = find(Ud(bmu_geo,1:munits)>kradius);

    Dst=dist( D(i,3:end),M(:,3:end)');
    Dst(non_searchable_units)=inf;
    [qerror_dta(i) bmu_dta]=min(Dst);
    
    Bmus(i)=bmu_dta;
end;
