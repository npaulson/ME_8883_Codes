
function [sphxo, A, dim] = statinfo(blen,ra,box,chids,bnds,loc,tme)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

color = hsv(20);

% number of large voxels in x,y and z
dim(1) = floor((2*bnds(1,2))/blen(1));
dim(2) = floor((2*bnds(2,2))/blen(2));
dim(3) = floor((2*bnds(3,2))/blen(3));

val = 1:length(loc(:,1));

chain_sep = accumarray(chids,val,[],@(x) {x});

for chid = 1:20

    chain = loc(chain_sep{chid},:);
    
    % 3d plot of all chains
    figure(1)
    
    H = plot3(chain(:,3), chain(:,1), chain(:,2));
    set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
        'MarkerFaceColor',color(chid,:),'MarkerSize',5)
    axis tight equal; grid on;    
    title2 = ['All Chains, Timestep = ', int2str(tme)];
    title(title2)

    if chid == 1; hold on; end;
    
end

subs(:,1) = floor(loc(:,1) / blen(1)) + 1;
subs(:,2) = floor(loc(:,2) / blen(2)) + 1;
subs(:,3) = floor(loc(:,3) / blen(3)) + 1;

A = accumarray(subs,val,[],@(x) {x});
A(end,:,:) = [];
A(:,end,:) = [];
A(:,:,end) = [];

cornercells = [   1,   1,   1;...
                  1,   1,dim(3);...
                  1,dim(2),   1;...
                  1,dim(2),dim(3);...
               dim(1),   1,   1;...
               dim(1),dim(2),   1;...
               dim(1),   1,dim(3);...
               dim(1),dim(2),dim(3)];

for mm = 1:8
    binloc = cornercells(mm,:) - 1;
    xs = blen(1)*(box(:,1) + binloc(1) + 0.5);
    ys = blen(2)*(box(:,2) + binloc(2) + 0.5);
    zs = blen(3)*(box(:,3) + binloc(3) + 0.5);

    plot3(zs,xs,ys,...
        'k-','LineWidth',2.0)
    axis tight equal; grid on;  
end

% here we develop the shape of the monomer
linvec = linspace(1,2*ra+1,2*ra+1);
[X,Y,Z] = meshgrid(linvec,linvec,linvec);
distances = sqrt((X-(ra + 1)).^2 + (Y-(ra + 1)).^2 + (Z-(ra + 1)).^2);

sphindxl = find(distances <= (ra + 0.5));

sphxo = zeros(length(sphindxl),3);
[sphxo(:,1), sphxo(:,2), sphxo(:,3)] = ind2sub(size(distances), sphindxl);

% plot the shape of the monomer
figure(2)
scatter3(sphxo(:,1), sphxo(:,2), sphxo(:,3))
axis equal

end