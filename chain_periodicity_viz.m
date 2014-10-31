% This script computes and plots the pair correlation function for one
% chain versus another at the specified timestep. The two chains are also
% plotted in 3D. The first chain is the 'reference chain' and is plotted
% once, the second is repeated periodically. The pair correlation is only
% computed within a specified radius for each of the monomers making up the
% reference chain.
%%% Noah Paulson, 9/30/2014

if exist('timestep') == 0
%     load total_info_11.mat
    load tilt6l.mat
end

close all

color = hsv(20);

ti = 1;

%2 is good

ch1 = 2; % ID of first chain
ch2 = 2; % ID of second chain


chain1 = (locations(:,2,ti) == ch1);
x_chain1 = chain1 .* locations(:,4,ti);
y_chain1 = chain1 .* locations(:,5,ti);
z_chain1 = chain1 .* locations(:,6,ti);
x_chain1(x_chain1==0) = [];
y_chain1(y_chain1==0) = [];
z_chain1(z_chain1==0) = [];

chains_int = [x_chain1, y_chain1, z_chain1];
ch_len = length(chains_int(:));

chain2 = (locations(:,2,ti) == ch2);
x_chain2 = chain2 .* locations(:,4,ti);
y_chain2 = chain2 .* locations(:,5,ti);
z_chain2 = chain2 .* locations(:,6,ti);
x_chain2(x_chain2==0) = [];
y_chain2(y_chain2==0) = [];
z_chain2(z_chain2==0) = [];

chains_per = [];

% the following code generates a specified tiling of simulation
% volumes (according to the periodic boundary conditions)
for xx = [0,1]        
    for yy = [-1,0]
        for zz = [-1,0,1]
            
            if xx == 0 && yy == 0 && zz == 0 
            
            else
                
                x_chain_per = x_chain2 + xx*2*bounds(1,2,ti)*ones(size(x_chain2));
                y_chain_per = y_chain2 + yy*2*bounds(2,2,ti)*ones(size(y_chain2));
                z_chain_per = z_chain2 + zz*2*bounds(3,2,ti)*ones(size(z_chain2));

                chains_spec = [x_chain_per,y_chain_per,z_chain_per];    
                chains_per = [chains_per ; chains_spec];
            
            end
        end
    end
end

% 3d plot of the specific chain of interest
figure(1)

clf

H = plot3(chains_int(:,3), chains_int(:,1), chains_int(:,2));
set(H,'LineStyle','none','Marker','o','MarkerEdgeColor',0.2*color(19,:),...
    'MarkerFaceColor',color(19,:),'MarkerSize',5)
hold on
% J = plot3(chains_per(:,1), chains_per(:,2), chains_per(:,3));
% set(J,'LineStyle','none','Marker','o','MarkerEdgeColor',.1*[1,1,1],...
%     'MarkerFaceColor','none','MarkerSize',5)
J = plot3(chains_per(:,3), chains_per(:,1), chains_per(:,2));
set(J,'LineStyle','none','Marker','o','MarkerEdgeColor',0.3*color(8,:),...
    'MarkerFaceColor',color(8,:),'MarkerSize',5)

sql = [bounds(1,1,ti),bounds(2,1,ti),bounds(3,1,ti);...
    bounds(1,2,ti),bounds(2,1,ti),bounds(3,1,ti);...
    bounds(1,2,ti),bounds(2,2,ti),bounds(3,1,ti);...
    bounds(1,1,ti),bounds(2,2,ti),bounds(3,1,ti);...
    bounds(1,1,ti),bounds(2,1,ti),bounds(3,1,ti)];
squ = [bounds(1,1,ti),bounds(2,1,ti),bounds(3,2,ti);...
    bounds(1,2,ti),bounds(2,1,ti),bounds(3,2,ti);...
    bounds(1,2,ti),bounds(2,2,ti),bounds(3,2,ti);...
    bounds(1,1,ti),bounds(2,2,ti),bounds(3,2,ti);...
    bounds(1,1,ti),bounds(2,1,ti),bounds(3,2,ti)];
edge1 = [bounds(1,1,ti),bounds(2,1,ti),bounds(3,1,ti);...
    bounds(1,1,ti),bounds(2,1,ti),bounds(3,2,ti)];
edge2 = [bounds(1,2,ti),bounds(2,1,ti),bounds(3,1,ti);...
    bounds(1,2,ti),bounds(2,1,ti),bounds(3,2,ti)];
edge3 = [bounds(1,1,ti),bounds(2,2,ti),bounds(3,1,ti);...
    bounds(1,1,ti),bounds(2,2,ti),bounds(3,2,ti)];
edge4 = [bounds(1,2,ti),bounds(2,2,ti),bounds(3,1,ti);...
    bounds(1,2,ti),bounds(2,2,ti),bounds(3,2,ti)];

box = {sql,squ,edge1,edge2,edge3,edge4};

for m = 1:length(box)
    plot3(box{m}(:,3),box{m}(:,1),box{m}(:,2),...
        'b--','LineWidth',1.5)
%         'LineStyle','--','LineWidth',1.5,'LineColor','c')
end
hold off

axis tight equal
grid on;    
title2 = ['Plot of chain ', int2str(ch1), ', including periodicity, Timestep = ', num2str(timestep(ti))];
title(title2);



