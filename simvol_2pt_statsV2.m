if exist('timestep','var') == 0
    load tilt6l.mat
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


loc = locations(:,4:6,ti);
vloc = ones(size(loc(:,1)));
translation = [bounds(1,1,ti)*vloc,bounds(2,1,ti)*vloc,bounds(3,1,ti)*vloc];
loc = loc - translation;


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
    
%     pause(.1)
    if chid == 1; hold on; end;
    
end

scale = 2;
locB = scale*loc + 1;
locC = round(locB);

for vv = 1:3
    el(vv) = ceil(scale*2*bounds(vv,2,ti));
end

ms = zeros(el(1),el(2),el(3));
for aa = 1 : length(vloc)
    ms(locC(aa,1),locC(aa,2),locC(aa,3)) = 1;
    
%     figure(2)
%     % plot the box which represents the closest bin to each particle
%     if mod(aa,100) == 0
%         xs = box(:,1) + locC(aa,1);
%         ys = box(:,2) + locC(aa,2);
%         zs = box(:,3) + locC(aa,3);
%         plot3(xs,ys,zs,...
%             'b-','LineWidth',0.5)
%         axis equal; grid on;
%         axis([0.5 (el(1) + 0.5)...
%               0.5 (el(2) + 0.5)...
%               0.5 (el(3) + 0.5)])
%         hold on
%     end
end

% plot3(locC(:,1),locC(:,2),locC(:,3),...
%     'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
%     'MarkerFaceColor',color(10,:),'MarkerSize',5);
% axis equal; grid on;
% axis([0.5 (el(1) + 0.5)...
%       0.5 (el(2) + 0.5)...
%       0.5 (el(3) + 0.5)])
% hold off

% Calculate Autocorrelation of particles
[T, xx] = SpatialStatsFFT(ms,[],'display',false,'shift',true);
% [T, xx] = SpatialStatsFFT(ms,[],'display',false,'shift',true,'cutoff',20);


% Plot point where 2pt stats are non-zero
Timg = find(T >= 4.5E-5);
sc = zeros(length(Timg),3);
[sc(:,1), sc(:,2), sc(:,3)] = ind2sub(size(T), Timg);

% xpos = (sc(:,1) >= 0 & sc(:,1) <= el(1));
% ypos = (sc(:,2) >= 0 & sc(:,2) <= el(2));
% zpos = (sc(:,3) >= floor(0.5*el(3) - 10) & sc(:,3) <= ceil(0.5*el(3) + 10));
% mult0 = xpos .* ypos .* zpos;
% 
% xvox = mult0 .* sc(:,1);
% yvox = mult0 .* sc(:,2);
% zvox = mult0 .* sc(:,3);
% xvox(xvox==0) = [];
% yvox(yvox==0) = [];
% zvox(zvox==0) = []; 

figure(3)
% plot3(xvox,yvox,zvox,...
% subplot(1,2,1)
plot3(sc(:,1),sc(:,2),sc(:,3),...
    'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
    'MarkerFaceColor',color(17,:),'MarkerSize',4);
axis equal tight;
title('2-pt statistics, probabilities > 4.5E-5');

% subplot(1,2,2)
% plot3(sc(:,1),sc(:,2),sc(:,3),...
%     'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
%     'MarkerFaceColor',color(17,:),'MarkerSize',4);
% axis equal tight;

% % Plot slice of 2-pt stats as image
% figure(4)
% % image(T(:,ceil(0.5*el(1)),:),'CDataMapping','scaled')
% image(T(:,:,ceil(0.5*el(3)+1)),'CDataMapping','scaled')
% colormap('jet')
% axis equal tight;
% colorbar
% shading flat
% caxis([0 4E-4])
% title('2-pt statistics, middle z slice');


% plot of the color choices
figure(5)
for cc = 1:20
    plot(cc,5,...
        'LineStyle','none','Marker','o','MarkerEdgeColor','k',...
        'MarkerFaceColor',color(cc,:),'MarkerSize',7);
    hold on
end
hold off