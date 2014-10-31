for ii = 1:6    
    close(figure(ii))
end


filename = 'tilt6l.mat';

if exist('timestep','var') == 0
    load(filename)
end

load plotbox

ti = 1;
ra = 4;
% roughly number of voxels per angstrom
scale = 6;

bnds = bounds(:,:,ti);
chids = locations(:,2,ti);
tme = timestep(ti);

% coordinates of all monomers
loc = locations(:,4:6,ti);
vloc = ones(size(loc(:,1)));
translation = [bounds(1,1,ti)*vloc,bounds(2,1,ti)*vloc,bounds(3,1,ti)*vloc];
loc = loc - translation;

blen = [35,9,9];

[sphxo, A, dim] = statinfo(blen,ra,box,chids,bnds,loc,tme);

density = zeros(dim(2),dim(3));
TimgL = zeros(dim(2),dim(3));

tic

vvec = [1,3,5];
voxstuff(box, blen, bnds, sphxo, 1, loc, A, vvec, scale, ra)

for xx = 1 : dim(2)
    for yy = 1 : dim(3)
        vvec = [1,xx,yy];
        [ density(xx,yy), ms, T, Timg ] = voxstuff(box, blen, bnds, sphxo, 0, loc, A, vvec, scale, ra);
        TimgL(xx,yy) = length(Timg);
        disp([xx,yy])
    end
end

toc

% Plot density
figure(7)
image(density,'CDataMapping','scaled')
colormap('jet')
axis equal tight;
colorbar
shading flat
title('Density Map')

% Plot number of statistics exceeding specified threshold
figure(8)
image(TimgL./density,'CDataMapping','scaled')
colormap('jet')
axis equal tight;
colorbar
shading flat
title('Crystalline Probability Map')

