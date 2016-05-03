function polylines=splitPolylines(multipoly)
%SPLITPOLYLINES Splits NaN separated polylines into a cell vector of polylines.
% example: multipoly=[1 1;1 2;2 2;2 1;NaN NaN;1.5 1.5;1.5 2.5;2.5 2.5;2.5 1.5]
% polygons=split_polygons(multipoly)
% Author: Florian de Boissieu
% Created: 2016-05-02

if ndims(multipoly)~=2 || size(multipoly,2)~=2
    error('Wrong input dimension: ''multipoly'' must be a N*2 matrix');
end
if ~any(isnan(multipoly(:,1)))
  error('Less than 2 polygons in multipoly');
end
    
i_poly=[0; find(isnan(multipoly(:,1)))];

polygons=cell();

if length(i_poly)>1
  for i=2:length(i_poly)
    polygons{i-1}=multipoly(i_poly(i-1)+1:i_poly(i)-1,:);
  end
end


    