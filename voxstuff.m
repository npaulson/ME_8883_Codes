function [ density, ms, T, Timg ] = voxstuff(box, blen, bnds, sphxo, pltfig, loc, A, vvec, scale, ra)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

vox0 = loc(A{vvec(1),vvec(2),vvec(3)},:);
voxv = ones(size(vox0(:,1)));
translation = [(vvec(1)-1)*blen(1)*voxv,...
               (vvec(2)-1)*blen(2)*voxv,...
               (vvec(3)-1)*blen(3)*voxv];
vox0b = vox0 - translation;
vox0c = scale*vox0b + 1;
vox0d = round(vox0c);
elx = blen(1)*scale + 1;
ely = blen(2)*scale + 1;
elz = blen(3)*scale + 1;

% figure(190)
% plot3(vox0(:,3), vox0(:,1), vox0(:,2),...
% 'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
% 'MarkerFaceColor',color(17,:),'MarkerSize',5); 
% axis equal
% axis([0 2*bnds(3,2) ...
%       0 2*bnds(1,2) ...
%       0 2*bnds(2,2)])

ms = zeros(elx+2*ra,ely+2*ra,elz+2*ra);

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

% Calculate Autocorrelation of particles
% [T, xx] = SpatialStatsFFT(ms,[],'display',false,'cutoff',blen(3),'shift',true);
[T, ~] = SpatialStatsFFT(ms,[],'display',false,'shift',true);

Timg = find(T > 0.008);

density = sum(ms(:))/length(ms(:));
    
if pltfig == 1

    color = hsv(20);
    
    % Plot box around the voxel of interest in the original image
    figure(1)

    load plotbox
    xs = blen(1)*(box(:,1) + vvec(1) - 0.5);
    ys = blen(2)*(box(:,2) + vvec(2) - 0.5);
    zs = blen(3)*(box(:,3) + vvec(3) - 0.5);

    hold on
    plot3(zs,xs,ys,...
        'c-','LineWidth',2.0)
    axis tight equal; grid on; 
    hold off


    figure(3)

    mspt = find(ms == 1);
    msplt = zeros(length(mspt),3);
    [msplt(:,1), msplt(:,2), msplt(:,3)] = ind2sub(size(ms), mspt);
    plot3(msplt(:,3), msplt(:,1), msplt(:,2),...
            'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
            'MarkerFaceColor',color(17,:),'MarkerSize',4); 
    hold on
    plot3(vox0c(:,3),vox0c(:,1),vox0c(:,2),...
        'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
        'MarkerFaceColor',color(5,:),'MarkerSize',10);

    axis equal; grid on;
    axis([0.5 (blen(3)*scale + 0.5)...
          0.5 (blen(1)*scale + 0.5)...
          0.5 (blen(2)*scale + 0.5)])
    title(['Monomers within voxel # ', int2str(sub2ind(size(A),vvec(1),vvec(2),vvec(3)))])
    hold off


    % Plot points where 2pt stats are greater than a certain value
    figure(4)

    sc = zeros(length(Timg),3);
    [sc(:,1), sc(:,2), sc(:,3)] = ind2sub(size(T), Timg);

    plot3(sc(:,1),sc(:,2),sc(:,3),...
        'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
        'MarkerFaceColor',color(17,:),'MarkerSize',5);
    axis tight equal;
    title('2-point statistics for individual voxel')


    % Plot slice of 2-pt stats as image
    figure(5)
    image(T(:,:,ceil(0.5*elz)),'CDataMapping','scaled')
    % image(squeeze(T(ceil(0.5*elx),:,:)),'CDataMapping','scaled')
    colormap('jet')
    axis equal tight;
    colorbar
    shading flat
    title('Slice of 2-pt statistics')


    % Plot histogram of 2-pt statistics probabilities
    figure(6)
    intv = 0.02;
    binranges = 0:intv:1;
    bincenters = binranges + 0.5*intv;
    hc = histc(T(:)/density,binranges);
    plot(bincenters,hc)
    axis([0 1 0 3000])
    xlabel('Normalized Autocorrelation Probability'); ylabel('Normalized Bincount');
    title('Histogram of Autocorrelation Probabilities')
end

end