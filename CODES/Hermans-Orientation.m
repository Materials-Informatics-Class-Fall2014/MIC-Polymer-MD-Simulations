clear;
%This program will plot a line between each successive
%atom in the chain. It then takes those vectors and
%gets the angle between it and the z-axis, for the 
%Herman's orientation.

openfile = '..\Polymer-MD\Converted\';
pictures = '..\Polymer-MD\Images\';
file = 'dump.deform375_10e7_5aniso_14';

D = dir([openfile, '*.txt']);
NumFiles = length(D(not([D.isdir])));
l = 1;
fid = fopen(horzcat(file,'_300.txt'));
fgetl(fid);
time = fscanf(fid, '%d\n');

for ii = 1:3
    fgetl(fid);
end

boundX = fscanf(fid, '%f %f\n', [1 2]);
boundY = fscanf(fid, '%f %f\n', [1 2]);
boundZ = fscanf(fid, '%f %f\n', [1 2]);
xvoxel = linspace(boundX(1),boundX(2), 2);
yvoxel = linspace(boundY(1),boundY(2), 2);
zvoxel = linspace(boundZ(1),boundZ(2), 13);

divisionsX = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
divisionsY = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
divisionsZ = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
lines = {0,0,0};
fgetl(fid);

data = fscanf(fid,'%d %d %*d %f %f %f\n', [5 7999])';
k = 1;
while k <= numel(data(:,2))
    j = 1;
    for j = 1:max(data(:,2))
        if data(k,2) == j
            divisionsX{j} = horzcat(divisionsX{j},[data(k,3)]);
            divisionsY{j} = horzcat(divisionsY{j},[data(k,4)]);
            divisionsZ{j} = horzcat(divisionsZ{j},[data(k,5)]);
            break;
        else
            j = j+1;
        end  
    end
    k = k+1;
end


H = plot3(divisionsX{1},divisionsY{1},divisionsZ{1});
set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',12);
axis tight equal;
grid on;
hold on;
S = plot3(divisionsX{2},divisionsY{2},divisionsZ{2});
set(S,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',12);
%Parse data according to voxel and store all data for each voxel in new
%array
%After that go through each voxels data and calculate the best fit lines
%and store the one best fit line through all of them.
numcells = (numel(xvoxel)-1)*(numel(yvoxel)-1)*(numel(zvoxel)-1);
ParsedData = cell(numcells,1);

for i = 1:numel(divisionsX)
    
    x = divisionsX{i}';
    y = divisionsY{i}';
    z = divisionsZ{i}';
    
    for h = 1:numel(x)
        p = 1;
        broken = 0;
        for j = 1:(numel(xvoxel)-1)
            for k = 1:(numel(yvoxel)-1)
                for l = 1:(numel(zvoxel)-1)
                    if (x(h)>xvoxel(j) && x(h)<xvoxel(j+1)&&...
                            y(h)>yvoxel(k) && y(h)<yvoxel(k+1)&&...
                            z(h)>zvoxel(l) && z(h)<zvoxel(l+1))
                        ParsedData{p} = vertcat(ParsedData{p},[i,x(h),y(h),z(h)]);
                        broken = 1;
                        break;
                    else
                        p = p+1;
                    end
                end
                if broken == 1
                    break;
                end
            end
            if broken == 1
                break
            end 
        end
    end
end

[n,mx] = size(divisionsX{1}); [ny,my] = size(divisionsY{1}); [nz,mz] = size(divisionsZ{1});
m = [mean(divisionsX{1}), mean(divisionsY{1}), mean(divisionsZ{1})];
w = [divisionsX{1}'-m(1), divisionsY{1}'-m(2), divisionsZ{1}'-m(3)];
a = (1/n)*w'*w;
[u,d,v] = svd(a);
q = u(:,1);
s = d(2,2)+d(3,3);

t = linspace(-100,100);
x = m(1)+q(1)*t;
y = m(2)+q(2)*t;
z = m(3)+q(3)*t;
plot3(x,y,z, 'LineWidth',5);

[n,mx] = size(divisionsX{2}); [ny,my] = size(divisionsY{2}); [nz,mz] = size(divisionsZ{2});
m = [mean(divisionsX{2}), mean(divisionsY{2}), mean(divisionsZ{2})];
w = [divisionsX{2}'-m(1), divisionsY{2}'-m(2), divisionsZ{2}'-m(3)];
a = (1/n)*w'*w;
[u,d,v] = svd(a);
p = u(:,1);
s = d(2,2)+d(3,3);

angle = atan2(norm(cross(q,p)),dot(q,p));
Hermans = (cos(angle)).^2
x1 = m(1)+q(1)*t;
y1 = m(2)+q(2)*t;
z1 = m(3)+q(3)*t;
plot3(x1,y1,z1, 'LineWidth',5);

xlim([boundX(1) boundX(2)]);
ylim([boundY(1) boundY(2)]);
zlim([boundZ(1) boundZ(2)]);

%Now we want to iterate through ParsedData, at each voxel calculate line of
%best fit for each separate chain, then line of best fit through all of
%them.
startpoint = cell(numcells,1);
direction = cell(numcells,1);
for i = 1:numel(ParsedData)
    %Iterate through voxel and extract x,y,z info based on chain id
    x = cell(20,1);
    y = cell(20,1);
    z = cell(20,1);
    
    for j = 1:(numel(ParsedData{i})/4)
        for k = 1:20
            if k == ParsedData{i}(j,1)
                %divide data in each voxel into chain id
                x{k} = vertcat(x{k},ParsedData{i}(j,2));
                y{k} = vertcat(y{k},ParsedData{i}(j,3));
                z{k} = vertcat(z{k},ParsedData{i}(j,4));
            end
        end
    end
    
    for l = 1:numel(x)
        if isempty(x{l}) == 0
            [n,mx] = size(x{l}); [ny,my] = size(y{l}); [nz,mz] = size(z{l});
            m = [mean(x{l}), mean(y{l}), mean(z{l})];
            w = [x{l}-m(1), y{l}-m(2), z{l}-m(3)];
            a = (1/n)*w'*w;
            [u,d,v] = svd(a);
            p = u(:,1);
            s = d(2,2)+d(3,3);

            startpoint{i} = vertcat(startpoint{i},m);
            direction{i} = vertcat(direction{i},p');
        end
    end
end
%We now have, for each voxel, the starting points and directions of the
%chains that are contained within that voxel.
%We can now iterate through each voxel to get the average of the best fit
%line and then find the angle for each of the previous best fit lines for
%each chain.
AllAngles = cell(numcells,1);
for i = 1:numel(direction)
    x = sum(direction{i}(:,1))/numel(direction{i}(:,1));
    y = sum(direction{i}(:,2))/numel(direction{i}(:,2));
    z = sum(direction{i}(:,3))/numel(direction{i}(:,3));
    xstart = sum(startpoint{i}(:,1))/numel(startpoint{i}(:,1));
    ystart = sum(startpoint{i}(:,2))/numel(startpoint{i}(:,2));
    zstart = sum(startpoint{i}(:,3))/numel(startpoint{i}(:,3));
    p = [x,y,z]./norm([x,y,z]);
    m = [xstart,ystart,zstart];
    avgstart{i} = m;
    avgdir{i} = p;
    
    for j = 1:numel(direction{i}(:,1))
        angle = atan2(norm(cross(direction{i}(j,:),avgdir{i})),dot(direction{i}(j,:),avgdir{i}));
        %if angle > (pi/2)
         %   angle = pi - angle;
        %end
        AllAngles{i} = horzcat(AllAngles{i},angle);
    end
end

x = cell(numcells,1);
y = cell(numcells,1);
z = cell(numcells,1);

for i = 1:numel(avgdir)
    t = linspace(0,100);
    x{i} = avgstart{i}(1)+t.*avgdir{i}(1);
    y{i} = avgstart{i}(2)+t.*avgdir{i}(2);
    z{i} = avgstart{i}(3)+t.*avgdir{i}(3);
end

for h = 1:numel(x)
    for j = 1:(numel(xvoxel)-1)
        for k = 1:(numel(yvoxel)-1)
            for l = 1:(numel(zvoxel)-1)
                for m = 1:numel(x{h})
                    while (x{h}(m)>xvoxel(j) && x{h}(m)<xvoxel(j+1)&&...
                            y{h}(m)>yvoxel(k) && y{h}(m)<yvoxel(k+1)&&...
                            z{h}(m)>zvoxel(l) && z{h}(m)<zvoxel(l+1)) == true
                    
                        x{h}(m) = [];
                        y{h}(m) = [];
                        z{h}(m) = [];
                    end
                end
            end
        end
    end
end

for i = 1:numel(x)
    plot3(x{i}, y{i}, z{i});
    hold on;
    axis tight equal;
    grid on;
end

%We now have an array of the avg start and avg direction vector for all
%points in each voxel. And the angle of all chain best fits to that line of
%best fit. We can now go through and find the hermans orientation with
%error as well.

for i = 1:numel(AllAngles)
    
    for j = 1:numel(AllAngles{i})
        Herm(j) = (cos(AllAngles{i}(j))).^2;
    end
    
    Hermans(i) = mean(Herm);
    HermansStdDev(i) = std(Herm);
end


fclose(fid);
