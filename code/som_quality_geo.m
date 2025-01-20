function [mqe_geo,mqe_dta,uqe_geo,uqe_dta] = som_quality_geo(sMap, D, bmu_ind,want_q_err)
% SOM_QUALITY_GEO Calculate the mean quantization for GEO problems
%
% [mqe_geo,mqe_dta,uqe_geo,uqe_dta] = som_quality_geo(sMap, D, bmu_ind,want_q_err)
%
%  qe = som_quality(sMap,D);
%
%  INPUTS  
%   sMap     (struct) a map struct
%   D                 the data
%            (struct) a data struct
%            (matrix) a data matrix, size dlen x dim
%   bmu_ind  (vector) a vector with the indexes of the bmus
%   want_q_error (bool) a boolean value indicating
%
% OUTPUTS
%   mqe_geo       (scalar) mean quantization error in geo-coords
%   mqe_dta       (scalar) mean quantization error in non-geo coords
%   uqe_geo       (scalar) quantization error in geo-coords for EACH UNIT
%   uqe_dta       (scalar) quantization error in non-geo coords for E.U. 

% Version 1.0   12/6/2004 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check arguments

% input arguments
if nargin < 2, error('Not enough input arguments.'); end

if nargin < 4, want_q_err =-1; end;
    
% data
if isstruct(D), D = D.data; end
[dlen dim] = size(D);

%define geo and non-geo indeces
ind_geo=[1 2];
ind_dta=[3:dim];
qe_geo=zeros(dlen,1);
qe_dta=qe_geo;

for i=1:dlen
    qe_geo(i)=dist(D(i,ind_geo),sMap.codebook(bmu_ind(i),ind_geo)');
    qe_dta(i)=dist(D(i,ind_dta),sMap.codebook(bmu_ind(i),ind_dta)');
end;

mqe_geo=mean(qe_geo);
mqe_dta=mean(qe_dta);

uqe_geo=[];
uqe_dta=[];
if want_q_err ~= -1
    [nunits x]=size(sMap.codebook);
    uqe_geo=zeros(nunits,1);
    uqe_dta=uqe_geo;
    for i=1:nunits
        ind = find( bmu_ind == i);
        uqe_geo(i)=mean(qe_geo(ind));
        uqe_dta(i)=mean(qe_dta(ind));
    end;
end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


