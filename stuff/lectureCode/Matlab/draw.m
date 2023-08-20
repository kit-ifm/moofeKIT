function draw(ed)

pprops = {'r','linewidth',1};
   
x=ed(:,1:2:size(ed,2))';
y=ed(:,2:2:size(ed,2))';  

plot(x,y,pprops{:})

axis equal

% xc=[x ; x(1,:)]; 
% yc=[y ; y(1,:)];
% 
% plot(xc,yc,pprops{:}) 

hold on
