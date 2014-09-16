% This code will read in the text output files from the LAMMPS molecular
% dynamics simulations.
% It then creates a new file of organized data in a textfile in a new
% directory.

original_dir = ('..\Polymer_MD\');
conv_dir = ('..\Polymer-MD\Converted\');

d = dir(original_dir);
file = 'dump.deform375_10e7_5aniso_11';
%for ii = find(~[dd.isdir])
 %   fo = fopen(horzcat(original_dir,dd(ii).name),'r');
 fo = fopen( file , 'r');
 t_ct = 1;
 while t_ct < 1250 && ~feof(fo)
     writefile = fopen(horzcat(conv_dir,file,'_',num2str(t_ct),'.txt'),'w');
     timestep = fgetl(fo);
     fprintf(writefile,'%s\n',timestep);
     t = fscanf(fo, '%i\n', 1);
     fprintf(writefile,'%i\n', t);
     
     for jj = 1:3
         headers = fgetl(fo);
         fprintf(writefile,'%s\n', headers);
     end
     
     lims = fscanf(fo, '%f %f\n',[ 2 3]);
     fprintf(writefile, '%f %f\n',lims);
     
     headers = fgetl(fo);
     fprintf(writefile,'%s\n', headers);
     
     A = fscanf(fo, '%d %d %d %f %f %f\n',[6 7999]);
     fprintf(writefile,'%d %d %d %f %f %f\n', A);
     t_ct = t_ct + 1;
 end
fclose(fo)
