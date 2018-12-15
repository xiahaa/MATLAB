function hand = multiple_boxplot(data,xlab,Mlab,colors)
% data is a cell matrix of MxL where in each element there is a array of N
% length. M is how many data for the same group, L, how many groups.
%
% Optional:
% xlab is a cell array of strings of length L with the names of each
% group
%
% Mlab is a cell array of strings of length M
%
% colors is a Mx4 matrix with normalized RGBA colors for each M.
% check that data is ok.
if ~iscell(data)
    error('Input data is not even a cell array!');
end
% Get sizes
M=size(data,2);
L=size(data,1);
if nargin>=4
    if size(colors,2)~=M
        error('Wrong amount of colors!');
    end
end
if nargin>=2
    if length(xlab)~=L
        error('Wrong amount of X labels given');
    end
end
% Calculate the positions of the boxes
hh = 0.5;
positions=1:hh:M*L*hh+1+hh*L;
positions(1:M+1:end)=[];
% Extract data and label it in the group correctly
x=[];
group1=[];
group2=[];
for ii=1:L
    for jj=1:M
        aux=data{ii,jj};
        x=vertcat(x,aux(:));
%         group=vertcat(group,ones(size(aux(:)))*jj+(ii-1)*M);
        group1=vertcat(group1,ones(size(aux(:)))*ii);
        group2=vertcat(group2,ones(size(aux(:)))*jj);
    end
end
group = [group1 group2];
% Plot it
hand = boxplot(x, group, 'positions', positions, 'FactorSeparator', 1:2);
% Set the Xlabels
aux=reshape(positions,M,[]);
labelpos = sum(aux,1)./M;
set(gca,'xtick',labelpos)
if nargin>=2
    set(gca,'xticklabel',xlab);
else
    idx=1:L;
    set(gca,'xticklabel',strsplit(num2str(idx),' '));
end
    
% Get some colors
if nargin>=4
    cmap=colors;
else
    cmap = hsv(M);
    cmap=vertcat(cmap,ones(1,M)*0.5);
end
color=repmat(cmap, 1, L);
% Apply colors
h = findobj(gca,'Tag','Box');
for jj=1:length(h)
   patch(get(h(jj),'XData'),get(h(jj),'YData'),color(1:3,jj)','FaceAlpha',color(4,jj));
end
if nargin>=3 && ~isempty(Mlab)
    legend(fliplr(Mlab),'Interpreter','latex','FontSize',8,'Location', 'northwest');
end

hpltall = get(gca,'children');   % the boxplot group
hbxplt = hpltall(end);
hall = get(hbxplt,'children');  % the individual components
hsepln = hall(end-2+1:end);     % the separator lines
delete(hsepln(1));

hlw = findall(gca,'tag','Lower Whisker');
huw = findall(gca,'tag','Upper Whisker');
set(hlw,'linestyle','-');
set(huw,'linestyle','-');
hout = findall(gca,'tag','Outliers');
% set(hout,'marker','.');
delete(hout);
end