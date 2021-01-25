function dendroOfX(T,Tlabels,n,font,X)
% replace T with a network as follows:
%   add a sink vertex such that each leaf has a directed edge to the sink
%   reverse the direction of all of the edges so the sink becomes the
%   source
commCoords=X(Tlabels,:);
%nodes=size(X,1);
%
x=X(1:n,1);y=X(1:n,2);
vertlabels = cellstr( num2str((1:n)') );  %' # labels correspond to their order

numclusts=size(Tlabels,2);

commlabels = cellstr( num2str((1:numclusts)') );

commx=commCoords(:,1);commy=commCoords(:,2);
% The code below creates the figure from the layout position data X
figure('units','normalized','position',[0.05 0.12 0.8 0.7]);
% Below we plot the edges colored by the edge weights stored in T
% first we make a suitable color map
mymap=[(0:1/127:1)',(1:-1/127:0)',ones(128,1)];
mymap=mymap.^6;
mymap(2:end,:)=mymap(2:end,:)/1.5;
wgPlotDendro(T,X,'edgeColorMap',mymap,'edgeWidth',1.2);
% axis square
hold on
plot(commx,commy,'s','MarkerEdgeColor','k','MarkerFaceColor',[1 0.8431 0], ...
    'MarkerSize',4);
text(x,y, vertlabels, 'VerticalAlignment','top', ...
                             'HorizontalAlignment','right','FontSize',font/1.5,...
                             'Rotation',45);
text(commx,commy, commlabels, 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right','FontSize',font,...
                             'Rotation',45);                             
% below we do not necessarily want
% clear x y labels
minT=min(T(T>0));
maxT=max(T(T>0));
if ~isequal(minT,maxT)
    step=(maxT-minT)/39;
    ticks=minT:step:maxT;
    colormap(mymap);
    caxis([minT maxT]);
    colorbar('Ticks',ticks);
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                BELOW ARE CALLED BY THE ABOVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hE]=wgPlotDendro(adjMat,coord,varargin)
%function [hE,hV]=wgPlot(adjMat,coord,varargin)
%
% Weighted Graph Plot from adjacency matrix [adjMat] and vertices
% coordinate [coord].
%
% INPUT:
%    [adjMat] = N x N sparse square adjacency matrix. 
%     [coord] = N x 2 matrix of vertex coordinates of the graph to be plotted.
%  [varargin] = are specified as parameter value pairs where the parameters are
%     'edgeColorMap' = m1 x 3 edge colormap for coloring edges by their 
%                      weight, {default=cool}
%        'edgeWidth' = scalar specifying the width of the edges. {default=0.1}
%     'vertexMarker' = single char, can be {.ox+*sdv^<>ph} as in plot {default='.'}
%   'vertexColorMap' = m2 x 3 vertex color map for coloring vertices by their
%                      weight, {default=summer}
%     'vertexWeight' = N x 1 vector of vertex weight that overrides using the
%                      diagonal of [adjMat] for specifying vertex size.
%      'vertexScale' = scalar vertex scaling factor for specifying scaling
%                      the size of vertices. {default=100}
%   'vertexmetadata' = N x 1 vector of vertex meta data for coloring the verticies.
% OUTPUT:
%  [hE] = vector/scalar of handles to edges of the drawn graph (1 per color).
%  [hV] = scalar handles to vertices of the drawn graph.
% 
% SEE ALSO: gplot, treeplot, spy, plot
%
% By: Michael Wu  --  michael.wu@lithium.com (May 2009)
%
%====================


% Set default parameter values
%--------------------
h=gca; prhDendro; hold on; axis off;
axesAreaDendro(h,[6 7 5 5]);
plotParm={'markerSize',5,'lineWidth',0.1,'marker','none','MarkerEdgeColor',[1,0.5,0.2]};
siz=size(adjMat);
vrtxSiz=100;
edgeMap=cool;
vrtxMap=summer;


% Parse parameter value pairs
%--------------------
nVarArgin=length(varargin);
for kk=1:2:nVarArgin
	switch lower(varargin{kk})
    case 'edgecolormap'
      edgeMap=varargin{kk+1};
    case 'edgewidth'
      plotParm=[plotParm,{'lineWidth',varargin{kk+1}}];
    case 'vertexmarker'
      plotParm=[plotParm,{'marker',varargin{kk+1}}];
    case 'vertexcolormap'
      vrtxMap=varargin{kk+1};
    case 'vertexweight'
      vrtxWt=varargin{kk+1};
    case 'vertexmetadata'
      vrtxCol=varargin{kk+1};
    case 'vertexscale'
      vrtxSiz=varargin{kk+1};
    otherwise
      error(['wgPlot >< Unknown parameter ''',varargin{kk},'''.']) ;
  end
end


% Determine if diagonal is weighted.
%--------------------
if exist('vrtxWt','var')
  vWt=vrtxWt;
else
  vWt=diag(adjMat);
end
vWeighted=length(setdiff(unique(vWt),0))>1;


% Map vWt weight to ?
%--------------------
if ~all(vWt==0)
  adjMat(speye(siz)~=0)=0;
end  % if ~any


% Determine if edges are weighted
%--------------------
[ii,jj,eWt] = find(adjMat);
qq=unique([ii,jj]);
minEWt=min(eWt);
maxEWt=max(eWt);
eWtRange=maxEWt-minEWt;
eWeighted=eWtRange>0;


% Map edge weight to edge colormap
%--------------------
if eWeighted
  neColor=size(edgeMap,1);
  eWt=ceil((neColor-1)*(eWt-minEWt)/(maxEWt-minEWt)+1);
end  % if eWtRange


% Plot edges
%--------------------
if eWeighted
  hE=[];
    
  for kk=1:neColor
    p=find(eWt==kk);
    nSegment=length(p);
    x=[coord(ii(p),1),coord(jj(p),1),repmat(nan,nSegment,1)]';
    y=[coord(ii(p),2),coord(jj(p),2),repmat(nan,nSegment,1)]';
    hE=[hE,plot(x(:),y(:),'color',edgeMap(kk,:),plotParm{:})];
    
  end  % for kk
  
 
else
    nSegment=length(ii);
    x=[coord(ii,1),coord(jj,1),nan(nSegment,1)]';
    y=[coord(ii,2),coord(jj,2),nan(nSegment,1)]';
    hE=plot(x(:),y(:),plotParm{:});
end  % if eWeighted


% Plot vertices
%--------------------



% Set axes
%--------------------
axis tight;
ax=axis;
dxRange=(ax(2)-ax(1))/500;
dyRange=(ax(4)-ax(3))/500;
axis([ax(1)-dxRange,ax(2)+30*dxRange,ax(3)-dyRange,ax(4)+dyRange]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H=prhDendro(H)
%function H=prh(H)
%
% sets the current figure to print on full page in LANDSCAPE mode.
%
% By Michael Wu  --  waftingpetal@yahoo.com (Oct 2001)
%
% ====================

if ~exist('H','var');
  H = gcf;
end;

papMarg=0.1;
set(H,'PaperOrientation','landscape','PaperPosition',[0+papMarg, 0+papMarg, 11-2*papMarg, 8.5-2*papMarg]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hAx]=axesAreaDendro(varargin)
%function [hAx]=axesArea(hAx,varargin)
%
% Set the margine of the axis by specifying a vector of distance to the figure edges.
%
% Input:
%  [varargin=hAx] = axis handle
%    [varargin=p] = position spec: 1, 2 or 4 element vector that specify the distance 
%                   from the edges of the figure by a percentage number between (0-49). 
%                   If 1 element it is used for all margins. 
%                   If 2 elements, p=[x-margins, y-margins].
%                   If 4 elements, p=[left, lower, right, upper] margins.
% Output:
%  [hAx] = axis handle of the axes.
%
% See also: axis, axes, ishandle, set
%
% By: Michael Wu  --  michael.wu@lithium.com (Mar 2009)
%
%====================


% Check if 1st input is axis handle
%--------------------
if ishghandle(varargin{1},'axes')
	hAx=varargin{1};
	p=varargin{2};
else
	hAx=gca;
	p=varargin{1};
end


% Process input arguments
%--------------------
p(p<0)=0;
p(p>49)=49;
p=p/100;


% Compute position property to be set
%--------------------
switch length(p)
	case 1
		xmin=p;
		ymin=xmin;
		xlen=1-2*p;
		ylen=xlen;
	case 2
		xmin=p(1);
		ymin=p(2);
		xlen=1-2*p(1);
		ylen=1-2*p(2);
	case 4
		xmin=p(1);
		ymin=p(2);
		xlen=1-p(1)-p(3);
		ylen=1-p(2)-p(4);	
	otherwise
		% Default Matlab position setting
		%--------------------
		xmin=0.13;
		ymin=0.11;
		xlen=0.775;
		ylen=0.815;	
end


% Set new position property of the axes
%--------------------
set(hAx,'position',[xmin ymin xlen ylen]);



