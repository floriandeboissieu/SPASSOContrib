function [firstX,indX,segX] = firstInterX(L1,L2)
%FIRSTINTERX first intersections of L2 over L1
%   [firstX, indX, segX] = firstInterX(L1,L2) returns the first intersection point of each L2 curve (or polyline) with any L1 
%   curve. The curves L1,L2 can be either closed or open and are described
%   by two-column-matrices, where each column contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]'
%
%   firstX is an Nx2 matrix, with N the number of curves of L2. Its columns correspond to the
%   x- and y- coordinates of the first intersection point of each L2 curve over any L1 curve. For each L2 curve not
%   intersecting any L1 curve, firstX row is set to NaN.
%   
%   IndX and SegX are Nx2 matrix giving respectivaly the index of the crossing segments within the multipolyline and within the polyline. 
%   First column is L1 index, second is L2 index. If L1 and L2 are only one curve each, SegX=IndX.
%   
%   Example:
%       L1=[1 1; 1 2; 2 2; 2 1; 1 1]; % square
%       L2=[0 0; 1.5 1.5; 0 3];
%       [firstX,indX,segX] = firstInterX(L1,L2);
%       
% Author: Florian de Boissieu
% Created: 2016-05-02,    using Octave 3.4.3
% Part of the code comes from InterX.m available on github/durka (credits : NS)

    hF = @le;
       
    %...Preliminary stuff
    x1  = L1(:,1);  x2 = L2(:,1)';
    y1  = L1(:,2);  y2 = L2(:,2)';
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    D1=D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1);
    D2=D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2');
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 

    %...Find multi-polylines nans
    inan=[0,find(isnan(x1) | isnan(y1))',Inf];
    jnan=[0,find(isnan(x2) | isnan(y2)),Inf];

    if isempty(i)
      firstX= nan(length(jnan)-1,2);
      indX= nan(length(jnan)-1,2);
      segX= nan(length(jnan)-1,2);
      return; 
    end
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = ([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L]);
    
    index=[i',j];
    jfirst=nan(length(jnan)-1,1);
    segX=nan(length(jnan)-1,2);
    %...Keep first of each 
    for jj=2:length(jnan)
      icurve=find(j>jnan(jj-1) & j<jnan(jj));
      if ~isempty(icurve)
        %jfirst(jj-1)=NaN;
      %else
        % find index of first crossing the point
        i1seg=find(j(icurve)==min(j(icurve))); % first crossing segment of icurve of L2
        [~,i1X]=min(sum((P(icurve(i1seg),:)-L2(j(icurve(i1seg)),:)).^2,2)); % first crossing point of crossing segment of L2
        jfirst(jj-1)=icurve(i1seg(i1X)); % points index in P
        segX(jj-1,2)=j(jfirst(jj-1))-jnan(jj-1);
        ii=find(i(jfirst(jj-1))>inan(1:end-1) & i(jfirst(jj-1))<inan(2:end));
        segX(jj-1,1)=i(jfirst(jj-1))-inan(ii);
      end
    end
    
    firstX=nan(length(jnan)-1,2);
    indX=nan(length(jnan)-1,2);
    firstX(~isnan(jfirst),:)=P(jfirst(~isnan(jfirst)),:); % change to Nx2 coordinate matrix
    indX(~isnan(jfirst),:)=index(jfirst(~isnan(jfirst)),:);
    
    
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end