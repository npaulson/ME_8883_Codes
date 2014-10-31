% This script computes and plots the pair correlation function for a single
% chain in a step of the MD simulation. It also plots the chain in 3D space
%%% Noah Paulson, 9/30/2014

if exist('timestep') == 0
    load tilt6l.mat
end

close all

color = hsv(20);

% ii: timestep id 
ii = 1;

c = 1;

for chid = 1:20

    chain = (locations(:,2,ii) == chid);
    x_chain = chain .* locations(:,4,ii);
    y_chain = chain .* locations(:,5,ii);
    z_chain = chain .* locations(:,6,ii);
    x_chain(x_chain==0) = [];
    y_chain(y_chain==0) = [];
    z_chain(z_chain==0) = [];

    chains_int = [x_chain, y_chain, z_chain];
    chains_ints{chid} = chains_int;
    ch_len = length(chains_int(:));
    
    chains_per = [];
    
    % the following code generates a specified tiling of simulation
    % volumes (according to the periodic boundary conditions)
    for xx = -1:1        
        for yy = -1:1
            for zz = -1:1
                x_chain_per = x_chain + xx*2*bounds(1,2,ii)*ones(size(x_chain));
                y_chain_per = y_chain + yy*2*bounds(2,2,ii)*ones(size(y_chain));
                z_chain_per = z_chain + zz*2*bounds(3,2,ii)*ones(size(z_chain));

                chains_spec = [x_chain_per,y_chain_per,z_chain_per];    
                chains_per = [chains_per ; chains_spec];
            end
        end
    end
    
    % range: maximum distance between particles to be considered
    range = 20;
    % dist: vector of the distances between all monomers
    [idx, dist] = rangesearch(chains_per,chains_int,range);
    
    dist_l = [];
    for rr = 1:length(dist)
        dist_l = [dist_l , cell2mat(dist(rr))];
    end
    dist_l = dist_l';
    dist_l(dist_l == 0) = [];
    
    intv = 0.2; % the interval in the binranges

    % binranges: vector of the starts/ends of the histogram bins
    binranges = 0 : intv : range;

    % bincenters: central radius value in each histogram bin
    bincenters = binranges(1:end) + 0.5 * intv;

    % bincounts: number of monomer distances per each bin
    bincounts(chid,:) = histc(dist_l,binranges)';
    
    % initialize pd_all_chains on first loop
    if c == 1
        bincounts_all_chains = zeros(1,length(bincounts(chid,:)));
    end
    
    bincounts_all_chains = bincounts_all_chains + bincounts(chid,:);
   
    c = c + 1;
end


colmean = mean(bincounts,1);

% if the mean of a column is exactly zero we assume all the numbers in the
% column are zero and set them to 1 to not screw up the normalization
colmean(colmean == 0) = 1;

colmean = repmat(colmean,length(bincounts(:,1)),1);

counts_colnorm = bincounts ./ colmean;

counts_colnorm(isnan(counts_colnorm)) = 1;

coeff = pca(counts_colnorm);

counts_pca = counts_colnorm * coeff;

pcadist = squareform(pdist(counts_pca));

compch = 17;
vecmax = [pcadist(compch,1:(compch-1)),0,pcadist(compch,(compch+1):end)];
maxindx = find(vecmax == max(vecmax));
vecmin = [pcadist(compch,1:(compch-1)),1000,pcadist(compch,(compch+1):end)];
minindx = find(vecmin == min(vecmin));

c = 1;

for chid = [compch, minindx, maxindx]
    
    figure(1)
    
    subplot(3,1,c)
    h2 = plot(binranges,bincounts_all_chains/(20), ...
    'b', 'LineSmoothing','on');
    hold on
    h1 = plot(binranges,bincounts(chid,:), ...
        'r', 'LineSmoothing','on');
    hold off

    axis([0 19 0 max(bincounts(:))])
    legend([h1,h2],'current chain','all chains')
    title1 = ['Pair correlation, chain ID = ', num2str(chid), ', timestep = ', num2str(timestep(ii))];
    title(title1)
    xlabel 'particle distance'
    ylabel 'pair count'

    figure(2)
    
    subplot(1,3,c)
    chainplt = chains_ints{chid};
    H = plot3(chainplt(:,1), chainplt(:,2), chainplt(:,3));
    set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(12,:),'MarkerSize',5)
    axis equal
    axis([bounds(1,1,ii) bounds(1,2,ii) bounds(2,1,ii) bounds(2,2,ii) bounds(3,1,ii) bounds(3,2,ii)])
    grid on;    
    title2 = ['Chain ',num2str(chid),', Timestep = ', num2str(timestep(ii))];
    title(title2)
    
    c = c + 1;
end

pcaA = 1;
pcaB = 2;
pcaC = 3;
labels = cellstr( num2str([1:20]') );

figure(3)

plot(counts_pca(:,pcaA),counts_pca(:,pcaB),'LineStyle','none','Marker','x')
text(counts_pca(:,pcaA),counts_pca(:,pcaB), labels, 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')
% text(counts_pca(:,pcaA),counts_pca(:,pcaB),counts_pca(:,pcaC), labels)

xlabel(['pca',int2str(pcaA)]); ylabel(['pca',int2str(pcaB)]); zlabel(['pca',int2str(pcaC)]);
title('Chains in PCA space via pair autocorrelations')
axis tight equal; grid on;

% plot3(counts_pca(:,pcaA),counts_pca(:,pcaB),counts_pca(:,pcaC),'LineStyle','none','Marker','x')
% % text(counts_pca(:,pcaA),counts_pca(:,pcaB),counts_pca(:,pcaC), labels, 'VerticalAlignment','bottom', ...
% %                              'HorizontalAlignment','right')
% text(counts_pca(:,pcaA),counts_pca(:,pcaB),counts_pca(:,pcaC), labels)
% 
% xlabel(['pca',int2str(pcaA)]); ylabel(['pca',int2str(pcaB)]); zlabel(['pca',int2str(pcaC)]);
% axis tight equal; grid on;
