clear
clc

if exist('timestep','var') == 0
    load tilt6l.mat
%     load perfect_xtal.mat
end

close all
color = hsv(20);
ti = 1;

% lines which plot a box centered at 0,0,0 with side lengths of 1 unit
box =   [ .5, .5,-.5;...
         -.5, .5,-.5;...
         -.5,-.5,-.5;...
          .5,-.5,-.5;...
          .5, .5,-.5;...
          .5, .5, .5;...
         -.5, .5, .5;...
         -.5,-.5, .5;...
          .5,-.5, .5;...
          .5, .5, .5;...
         -.5, .5, .5;...
         -.5, .5,-.5;...
         -.5, .5, .5;...
         -.5,-.5, .5;...
         -.5,-.5,-.5;... 
         -.5,-.5, .5;...
          .5,-.5, .5;...
          .5,-.5,-.5;...
          .5,-.5, .5]; 

% lengths of the large voxels in x, y and z
blenx = floor(2*bounds(1,2,ti));
bleny = floor(2*bounds(2,2,ti));
blenz = floor(2*bounds(3,2,ti));
% number of large voxels in x,y and z
xdim = floor((2*bounds(1,2,ti))/blenx);
ydim = floor((2*bounds(2,2,ti))/bleny);
zdim = floor((2*bounds(3,2,ti))/blenz);

% coordinates of all monomers
loc = locations(:,4:6,ti);

indc = zeros((xdim+1)*(ydim+1)*(zdim+1),3);

vloc = ones(size(loc(:,1)));
translation = [bounds(1,1,ti)*vloc,bounds(2,1,ti)*vloc,bounds(3,1,ti)*vloc];
loc = loc - translation;
% total number of large voxels
cmax =(xdim+1)*(ydim+1)*(zdim+1);


for chid = 1:20

    chain = (locations(:,2,ti) == chid);
    xvox = chain .* loc(:,1);
    yvox = chain .* loc(:,2);
    zvox = chain .* loc(:,3);
    xvox(xvox==0) = [];
    yvox(yvox==0) = [];
    zvox(zvox==0) = [];
    
    chains = [xvox,yvox,zvox];
    
    % 3d plot of all chains
    figure(1)
    
    H = plot3(chains(:,3), chains(:,1), chains(:,2));
    set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
        'MarkerFaceColor',color(chid,:),'MarkerSize',5)
    axis tight equal; grid on;    
    title2 = ['All Chains, Timestep = ', num2str(timestep(ti))];
    title(title2)

    if chid == 1; hold on; end;
    
end

tic
c = 1;
for xx = 0:xdim-1
    
    lx = xx * blenx;
    xpos = (loc(:,1) > lx & loc(:,1) < (lx + blenx));
    
    for yy = 0:ydim-1
        
        ly = yy * bleny;
        ypos = (loc(:,2) > ly & loc(:,2) < (ly + bleny));
        
        for zz = 0:zdim-1
        
            lz = zz * blenz;
            zpos = (loc(:,3) > lz & loc(:,3) < (lz + blenz));

            % vector of 1 and 0s where the one designates the preceding
            % conditions have been met. This will be multiplied by the
            % coordinates of the monomers to only get a monomer from  a
            % certain chain in a certain spatial bin.
            mult0 = xpos .* ypos .* zpos;

            xvox = mult0 .* loc(:,1);
            yvox = mult0 .* loc(:,2);
            zvox = mult0 .* loc(:,3);
            xvox(xvox==0) = [];
            yvox(yvox==0) = [];
            zvox(zvox==0) = []; 
            
            vox = [xvox,yvox,zvox];
            all_vox{c} = vox;
            
            indc(c,:) = [xx,yy,zz];
            c = c + 1;
           
            % Plot boxes around the 8 corners of the voxelized simulation
            % volume
            if (xx == 0 || xx == (xdim-1)) && ...
               (yy == 0 || yy == (ydim-1)) && ...
               (zz == 0 || zz == (zdim-1))
                
                xs = blenx*box(:,1) + lx + 0.5*blenx;
                ys = bleny*box(:,2) + ly + 0.5*bleny;
                zs = blenz*box(:,3) + lz + 0.5*blenz;
                
                plot3(zs,xs,ys,...
                    'k-','LineWidth',2.0)
                axis tight equal; grid on;     
                hold on
            end   
        end
    end
end
hold off

toc

% ID of the voxel of interest
vid = 1;
% roughly number of voxels per angstrom
scale = 3;

vox0 = all_vox{vid};
voxv = ones(size(vox0(:,1)));
translation = [indc(vid,1)*blenx*voxv,...
               indc(vid,2)*bleny*voxv,...
               indc(vid,3)*blenz*voxv];
vox0b = vox0 - translation;
vox0c = scale*vox0b + 1;
vox0d = round(vox0c);
elx = blenx*scale + 1;
ely = bleny*scale + 1;
elz = blenz*scale + 1;

% roughly the radius of the ethylene molecule
ra = 2;

ms = zeros(elx+2*ra,ely+2*ra,elz+2*ra);

% here we develop the shape of the monomer
linvec = linspace(1,2*ra+1,2*ra+1);
[X,Y,Z] = meshgrid(linvec,linvec,linvec);
distances = sqrt((X-(ra+1)).^2 + (Y-(ra+1)).^2 + (Z-(ra+1)).^2);

sphindxl = find(distances <= (ra + 0.5));

sphxo = zeros(length(sphindxl),3);
[sphxo(:,1), sphxo(:,2), sphxo(:,3)] = ind2sub(size(distances), sphindxl);

% plot the shape of the monomer
figure(2)
scatter3(sphxo(:,1), sphxo(:,2), sphxo(:,3))
axis equal



for aa = 1 : length(voxv)
    
    pos = vox0d(aa,:);
    sphx = sphxo;
    
    sphx(:,1) = sphx(:,1) - 1 + pos(1);
    sphx(:,2) = sphx(:,2) - 1 + pos(2);
    sphx(:,3) = sphx(:,3) - 1 + pos(3);

    sphxl = sub2ind(size(ms),sphx(:,1), sphx(:,2), sphx(:,3));
    ms(sphxl) = 1;
    
end

ms = ms(ra+1:end-ra, ...
        ra+1:end-ra, ...
        ra+1:end-ra);

% mspt = find(ms == 1);
% msplt = zeros(length(mspt),3);
% [msplt(:,1), msplt(:,2), msplt(:,3)] = ind2sub(size(ms), mspt);
% 
% figure(3)
% plot3(msplt(:,1), msplt(:,2), msplt(:,3),...
%         'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
%         'MarkerFaceColor',color(17,:),'MarkerSize',4); 
% hold on
% plot3(vox0c(:,1),vox0c(:,2),vox0c(:,3),...
%     'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
%     'MarkerFaceColor',color(5,:),'MarkerSize',10);
%     
% axis equal; grid on;
% axis([0.5 (blenx*scale + 0.5)...
%       0.5 (bleny*scale + 0.5)...
%       0.5 (blenz*scale + 0.5)])
% title(['Monomers within voxel # ', int2str(vid)])
% hold off


% Calculate Autocorrelation of particles
% [T, xx] = SpatialStatsFFT(ms,[],'display',false,'cutoff',blenz,'shift',true);
[T, xx] = SpatialStatsFFT(ms,[],'display',false,'shift',true);


% Plot points where 2pt stats are greater than a certain value
Timg = find(T > 0.025);
sc = zeros(length(Timg),3);
[sc(:,1), sc(:,2), sc(:,3)] = ind2sub(size(T), Timg);

figure(4)
plot3(sc(:,1),sc(:,2),sc(:,3),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
    'MarkerFaceColor',color(17,:),'MarkerSize',5);
axis tight equal;
title('2-point statistics for individual voxel')


% Plot slice of 2-pt stats as image
figure(5)
% image(T(:,:,ceil(0.5*elz)),'CDataMapping','scaled')
image(squeeze(T(ceil(0.5*elx),:,:)),'CDataMapping','scaled')
colormap('jet')
axis equal tight;
colorbar
shading flat
title('Slice of 2-pt statistics')
% caxis([0 1E-5])

% Plot histogram of 2-pt statistics probabilities
figure(6)
intv = 0.02;
binranges = 0:intv:1;
bincenters = binranges + 0.5*intv;
density = sum(ms(:))/length(ms(:));
hc = histc(T(:)/density,binranges);
plot(bincenters,hc)
axis([0.2 1 0 6000])
xlabel('Normalized Autocorrelation Probability'); ylabel('Normalized Bincount');
title('Histogram of Autocorrelation Probabilities')