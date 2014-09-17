clear;
openfile = '..\Polymer-MD\Converted\';
pictures = '..\Polymer-MD\Images\';
file = 'dump.deform375_10e7_5aniso_11';

%writerObj = VideoWriter('movie.avi');
%open(writerObj);

color = hsv(20);

D = dir([openfile, '*.txt']);
NumFiles = length(D(not([D.isdir])));
l = 1;

for jj = 1:50:NumFiles
    fid = fopen(horzcat(openfile,file,'_',num2str(jj),'.txt'), 'r');
    %fprintf('%d\n',jj)
    m = 1;
    while ~feof(fid)
        
        fgetl(fid);
        time = fscanf(fid, '%d\n');
        
        for ii = 1:3
            fgetl(fid);
        end
        
        boundX = fscanf(fid, '%f %f\n', [1 2]);
        boundY = fscanf(fid, '%f %f\n', [1 2]);
        boundZ = fscanf(fid, '%f %f\n', [1 2]);
        
        fgetl(fid);
        
        data = fscanf(fid,'%*d %d %*d %f %f %f\n', [4 7999])';
        
        k = 1;
        divisionsX = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        divisionsY = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        divisionsZ = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        
        while k <= numel(data(:,1))
            j = 1;
            for j = 1:max(data(:,1))
                if data(k,1) == j;
                    divisionsX{j} = horzcat(divisionsX{j},[data(k,2)]);
                    divisionsY{j} = horzcat(divisionsY{j},[data(k,3)]);
                    divisionsZ{j} = horzcat(divisionsZ{j},[data(k,4)]);
                    break;
                else
                    j = j+1;
                end
            end
            k = k + 1;
        end
            
            
        for m = 1:max(data(:,1))
            divisionsX{m}(1) = [];
            divisionsY{m}(1) = [];
            divisionsZ{m}(1) = [];
            
            
            H = plot3(divisionsX{m}, divisionsY{m}, divisionsZ{m});
       
            set(H,'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',color(m,:),'MarkerSize',12);
            %fprintf('%d %d %d\n', color(m,:))
            hold on;
            
        end
        
        %axis([min(data(:,2)), max(data(:,2)), min(data(:,3)), max(data(:,3)), min(data(:,4)), max(data(:,4))], 'equal')
        
        title(horzcat('Timestep = ', num2str(time)));
        grid on;
        hold off;
        axis tight equal
        saveas(H,horzcat(pictures,'dumpdeform375_10e7_5aniso_11','_',num2str(jj),'.png'));
        %xlim([-25 25]);
        %ylim([-15 15]);
        %zlim([-200 200]);
        %mov(l) = getframe;
        %writeVideo(writerObj,mov);
        l = l + 1;
        
    end
    fclose(fid);
   
end
 close(writerObj);
