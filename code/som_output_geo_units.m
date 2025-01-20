function outputData=som_output_geo_units( sM, sD, k) 
% function outputData=som_output_geo_units( sM, sD, k)
% Calculates a number of values concerning the neurons
%
% INPUTS
%   sM - Map for which we want to generate the output
%   sD - Data used to calculate the quantization errors of the units
%   k  - Geo_SOM radius used to calculate the Best Match (a value of -1
%   implies that a standard SOM will be used.
%
% OUTPUTS
%   outputData  Matrix with [n,5] elements. Each row is for a different unit.
%       column 1 - x coordinates
%       column 2 - y coordinates
%       column 3 - Geographical error for each unit
%       column 4 - Quantization error for each unit
%       column 5 - Average distance (in the data's NON-GEO dimension) to it's neighnors (u-MAT)
%
%   V.1.0 - By Lobo & Bacao, 2004/06/21

if k==-1                            % Procurar a BMU de acordo com as regras tradicioanis (k=-1)
    bmu=som_bmus(sM,sD);            % ou de acordo com as regras do GeoSOM
    k=inf;
else
    bmu=som_bmus_geo(sM,sD,k);
end;

[g d gm dm ]=som_quality_geo(sM,sD,bmu,k);

sU=som_umat_geo(sM);    % calculate UMAT-GEO
                        % Identify the values of the UMAT that correspond
                        % to neurons, and store them in a vector "sUlinear"

[maxx, maxy]=size(sU);
xi=1:2:maxx;yi=1:2:maxy;
sUunits=sU(xi,yi);      % use an intermediate 2D matrix with only the values
                        % of the UMAT that correspond to neurons

                        [maxx, maxy]=size(sUunits);
sUlinear=zeros(maxx*maxy,1);
for i = 1:(maxx*maxy),  % copy the values to the final (linear) vector
    sUlinear(i)=sUunits( floor((i-1)/maxy)+1,mod(i-1,maxy)+1 );
end;

outputData = [ sM.codebook(:,1) sM.codebook(:,2) gm dm sUlinear];
