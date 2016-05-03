function [touched,lons,lats,lonf,latf,touchedlon,touchedlat,touchedpoly]=timefrompolygon(day0,lons0,lats0,numdays,product,polygon,out2in)
% Computes the trajectory of the water mass observed on day0 at (lons0,lat0) until day0+numdays, 
% and extracts the time and coordinates of the first intersection with one of the polygon boundaries.
%
% Inputs:
%   day0: date vector. starting date for the advection.
%   lon0,lat0: coordinates of the watermass to advect.
%   numdays: number of backward days (if 'numdays'<0 it is a backward advection, otherwise it is made forward).
%   product: sea surface speed vectors (U,V). Possibilities are the ones in aviso_load.m from SPASSO package.
%   polygon: Nx2 matrix of coordinates of one a closed polyline. Multiple polygons can also be used separating them by a NaN row (ex:[1 1; 1 2; 2 2; 2 1; 1 1; NaN NaN; 3 4; 4 5; 1 5; 3 4])
%   out2in: 1 if only the crossing polygons from the outside is to take into account, 0 otherswise.
%
% Values:
%   touched: time that it took to touch a boundary, NaN if it did not touch.
%   lons, lats: an numdaysxN coordinates of the daily points of the backward advection derived from current spatial interpollations and Runge-Kutta algorithm.
%   lonf, latf: final coordiantes after numdays.
%   touchedlon,touchedlat: coordinates of the first intersection of the water mass trajectory with one of the polygons.
%   touchedpoly: index of the segments of polygon implied in the intersections.
%                polygon(touchedpoly,:) and polygon(touchedpoly+1,:) are the extremities of the intersecting segments.
%
% WARNING: numdays positive means a advection in future... should be homogenized with aviso_fsle() which is the oposit.
% WARNING: touched is 2-based instead of 1-based... should be modified with timefrombathy

if any(isnan(polygon)) % case of multi-polygon
  spolygon=splitPolylines(polygon); % split polygons, used for inpolygon function
else
  spolygon={polygon}; % case of a unique polygon
end

  
[lonf,latf,lons,lats]=aviso_lonlatadvect(day0,lons0,lats0,numdays,product); % advection of points from day0,lons0,lats0 to day0+numdays

lons=lons(1:8:end,:);% precision of 1 day, see aviso_lonlatadvect.m for details.
lats=lats(1:8:end,:);
slons=size(lons); % size of lons matrix
spoly=size(polygon);

maxproc=floor(2^21/(spoly(1)*slons(1))); % maximum number of lons,lats points that can be processed at the time, otherwise it leads to out of memory
nstep=ceil(slons(2)/maxproc); % number of computing steps to process all trajectories of the input points.


fprintf('\nComputing the intersections...');
if nstep==1
  % make a polyline from lons,lats to process at once in firstInterX.
  polyline=[reshape([lons;nan(1,size(lons,2))],[],1),reshape([lats;nan(1,size(lats,2))],[],1)]; % make it multipolyline Nx2
  polyline=polyline(1:end-1,:); % remove the trailing NaN
  [firstX,indX,segX]=firstInterX(polygon,polyline); % segX is the segment of each polyline which crosses polygon
else
  fprintf('\nTrajectory matrix too big, computing made in %d steps\n',nstep);
  firstX=nan(slons(2),2); % coordinates of the first intersection of each polyline of 'polyline' with any polygon of 'polygon'
  indX=nan(slons(2),2); % index of the segments implied in the intersection: first column is for 'polygon' variable, second column is for 'polyline' (i.e. lons,lats backward trajectory)
  segX=nan(slons(2),2); % 
  offset=0;
  for ii=1:nstep
      if ii==nstep
        istepll=(ii-1)*maxproc+1:slons(2);
      else
        istepll=(ii-1)*maxproc+1:ii*maxproc;
      end
      polyline=[reshape([lons(:,istepll);nan(1,length(istepll))],[],1),reshape([lats(:,istepll);nan(1,length(istepll))],[],1)]; % make it multipolyline Nx2
      polyline=polyline(1:end-1,:); % remove trailing NaN
      [firstX(istepll,:),indX(istepll,:),segX(istepll,:)]=firstInterX(polygon,polyline);
      indX(istepll,2)=indX(istepll,2)+offset;
      offset=size(polyline,1)+offset+1; % the +1 is for the NaN not counted
      % printf('Nstep=%03d\n',ii);
  end
  polyline=[reshape([lons;nan(1,size(lons,2))],[],1),reshape([lats;nan(1,size(lats,2))],[],1)]; % make it multipolyline Nx2
  polyline=polyline(1:end-1,:); % remove trailing NaN
end    
  
% detect the points that had an intersection with the polygons
plX=~isnan(indX(:,1));
% check if coming from inside a polygon (on the polygon is considered as inside)
inpoly=false(size(indX,1),1);
inpoly(plX)=any(cell2mat(cellfun(@(x) inpolygon(polyline(indX(plX,2),1),polyline(indX(plX,2),2),x(:,1),x(:,2)),spolygon,'UniformOutput',false)),2);

if out2in % set to NaN the intersections starting from inside the polygons
  firstX(inpoly,:)=NaN;
  indX(inpoly,:)=NaN;
  segX(inpoly,:)=NaN;
else % set to NanN the segments that are starting from outside or on the border of the polygons
    % check if coming from on a polygon (on the polygon is considered as inside)
  onpoly=false(length(indX),1);
  onpoly(plX)=any(cell2mat(cellfun(@(x) onpolygon(polyline(indX(plX,2),1),polyline(indX(plX,2),2),x(:,1),x(:,2)),spolygon,'UniformOutput',false)),2);
  firstX(~inpoly | onpoly,:)=NaN;
  indX(~inpoly | onpoly,:)=NaN;
  segX(~inpoly | onpoly,:)=NaN;
end

touchedlon=firstX(:,1);
touchedlat=firstX(:,2);
touched=segX(:,2);
touchedpoly=indX(:,1);

function on=onpolygon(X,Y,XV,YV)
% Proxy to inpolygon with 'on' only return
[~,on]=inpolygon(X,Y,XV,YV);

