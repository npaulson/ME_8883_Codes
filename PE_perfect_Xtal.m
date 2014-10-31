%This code develops perfect crystalline structure based on polyethylene
%information gathered from the Handbook of Polyethylene.

%For reference: PE forms a BC orthrhombic crystal structure

bond_angle = deg2rad(112);
bond_length = 1.53; %in angstroms

%initial dimensions of cell in simulation (angstroms)
xdim = 40;
ydim = 40;
zdim = 140;

%create a vector with the initial chain
%this vector can then be repeated at the unit cell corners to create a
%complete crystal over the volume
k = 1;
for ii = 0:1.27:zdim
    if rem(ii,2.54)==0
        initchain(k,1) = 0; %x location
        initchain(k,2) = 0; %y location
        initchain(k,3) = ii; %z location
        bc_chain(k,1) = 4.356;
        bc_chain(k,2) = 2.475;
        bc_chain(k,3) = ii;
    else
        initchain(k,1) = 0.646;
        initchain(k,2) = 0.561;
        initchain(k,3) = ii;
        bc_chain(k,1) = 3.71;
        bc_chain(k,2) = 3.037;
        bc_chain(k,3) = ii;
    end
    k=k+1;
end

%we now have one complete chain so we can repeat in a and b accordingly
%a rotated chain will need to be created for the (body-center)
chains = vertcat(initchain,bc_chain);
ctr = 1;
for ii = 7.42:7.42:xdim
    newchainx = initchain(:,1)+ctr*7.42;
    newbc_chainx = bc_chain(:,1)+ctr*7.42;
    chains = vertcat(chains,[newchainx,initchain(:,2),initchain(:,3)],[newbc_chainx,bc_chain(:,2),bc_chain(:,3)]);
    ctr = ctr + 1;
end;
ctr = 1;
allchains = chains;
for ii = 4.95:4.95:ydim
    newchainy = chains(:,2)+ctr*4.95;
    allchains = vertcat(allchains,[chains(:,1),newchainy,chains(:,3)]);
    ctr = ctr+1;
end
j = 1;
for ii = 1:numel(allchains(:,1))
    atoms(ii) = ii;
    id(ii) = j;
    if rem(ii,numel(initchain(:,1)))==0
        j = j+1; %will create unique chain id's
    end
    mol(ii) = 1;
end
num = numel(initchain(:,1));
col = numel(allchains(:,1))/num;
color = hsv(col);
H = plot3(allchains(1:num,1),allchains(1:num,2),allchains(1:num,3));
set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(1,:),'MarkerSize',5);
grid on;
hold on;
for m = 2:col
    H = plot3(allchains((m-1)*num+1:m*num,1),allchains((m-1)*num+1:m*num,2),allchains((m-1)*num+1:m*num,3));
    set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(m,:),'MarkerSize',5);
end
axis tight equal;
axis([0 xdim 0 ydim 0 zdim])

hold off;


% save in my format

locations = zeros(length(allchains),6,1);
locations(:,4:6,1) = allchains;

bounds = zeros(3,2,1);
bounds(1,:,1) = [0,xdim];
bounds(2,:,1) = [0,ydim];
bounds(3,:,1) = [0,zdim];

timestep = 1;

% save perfect_xtal.mat bounds locations timestep





