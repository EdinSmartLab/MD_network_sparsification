function [hE,hV]=weightedGraphPlot(adjMat,coord,vrtxWt)
% Plots the weighted graph inputs 1: Adjacency matrix, 2: vertex
% coordinates, 3: Vertex weights (for size of dot)

%% Initialize
hold on; axis off;
plotParm={'markerSize',6,'lineWidth',1.5,'marker','.','MarkerEdgeColor',[1,0.5,0.2]};
siz=size(adjMat);
vrtxSiz=5;
% edgeMap = flipud(parula);
% edgeMap = parula;
edgeMap=flipud(hot);
colormap(flipud(hot))%all
% edgeMap = flipud(copper);
% colormap(copper)

%% Vertex weights 
if exist('vrtxWt','var')
  vWt=vrtxWt;
else
  vWt=diag(adjMat);
end
vWeighted=length(setdiff(unique(vWt),0))>1;

%% Map edge weight to edge colormap
if ~all(vWt==0)
  adjMat(speye(siz)~=0)=0;
end  
[ii,jj,eWt] = find(adjMat);
qq=unique([ii,jj]);
minEWt=min(eWt);
maxEWt=max(eWt);
eWtRange=maxEWt-minEWt;
eWeighted=eWtRange>0;
if eWeighted
  neColor=size(edgeMap,1);
  eWt=ceil((neColor-1)*(eWt-minEWt)/(maxEWt-minEWt)+1);
end  
%% Plot edges
if eWeighted
  hE=[]; 
  neColor;
  for kk=1:neColor
    p=find(eWt==kk);
    nSegment=length(p);
    x=[coord(ii(p),1),coord(jj(p),1),repmat(nan,nSegment,1)]';
    y=[coord(ii(p),2),coord(jj(p),2),repmat(nan,nSegment,1)]';
    hE=[hE,plot(x(:),y(:),'color',edgeMap(kk,:),plotParm{:})];
  end  
else
    nSegment=length(ii);
    x=[coord(ii,1),coord(jj,1),nan(nSegment,1)]';
    y=[coord(ii,2),coord(jj,2),nan(nSegment,1)]';
    hE=plot(x(:),y(:),plotParm{:});
end 

%% Plot vertices
% if vWeighted
%   size(vWt);
%   minVwt=min(vWt);
%   maxVwt=max(vWt);
%   vWt=vrtxSiz*(vWt-minVwt)/(maxVwt-minVwt)+1;
%   hV=scatter(coord(qq,1),coord(qq,2),vWt,'filled','MarkerFaceColor',[1,0.5,0.2]);
% else
%   hV=plot(coord(qq,1),coord(qq,2),plotParm{:},'LineStyle','none');
% end
% colorbar;
% whitebg(1,'k');