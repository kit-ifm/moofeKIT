function draw_conf(Spann,q,edof,edofSpannung,ed)
%% Konturplot Spannung

[ed_neu]    = extract(edof,q);
[ed_Spann]  = extract(edofSpannung,Spann);

for i = 1:size(ed,1)
    X(1:4,i)    = ed_neu(i,[1 3 5 7])';
    Y(1:4,i)    = ed_neu(i,[2 4 6 8])';
    PK(1:4,i)   = ed_Spann(i,[1 2 3 4])';
end

patch(X,Y,PK)
caxis([0 max(Spann)]);
% shading interp
set(gcf,'color',[1 1 1])



