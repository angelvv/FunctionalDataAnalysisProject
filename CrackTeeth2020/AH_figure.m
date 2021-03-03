function fig = AH_figure(numRows, numCols, name)
%fig = AH_figure(numRows, numCols, name); %numRows, numCols, name

%screensize = get( groot, 'Screensize' );
if ~exist('name','var'); name = 'figure'; end
fig = figure('name',name,'Position',[5 30 320*numCols 270*numRows]);%bottom left x,y,width,height; y=0 can see title; y=30 can see x label
end
