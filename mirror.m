close all
clear all
 
%load data
path = '2_Reg_Sementations_to_1602'
directory =directory;
filesAndFolders = dir(directory);

case_name = "name";

save_directory = append(directory,'\',case_name,'\',case_name,'_','L_ROI_projected.stl');

Vessel_file_name = append(directory,'\',case_name,'\',case_name,'_Both.stl');
fardad1=stlread(Vessel_file_name);  %whole vessel

ROI_file_name = append(directory,'\',case_name,'\',case_name,'_L_ROI.stl');
fardad2=stlread(ROI_file_name);   %ROI

vrs = fardad1.Points;  %whole vessel
[m,n] = size(vrs);

ivrs = fardad2.Points; %isolated vessel

CM(1,1) = mean(vrs(:,1));
CM(2,1) = mean(vrs(:,2));
CM(3,1) = mean(vrs(:,3));

%median plane normal
CMM = repmat(CM',m,1);

centered_mean = vrs-CMM;

C = cov(centered_mean);

[EV,EC]=eig(C);   %PCA values (eigen values and eigen vectors of standardized matrix

NM = cross(EV(:,1),EV(:,2));

d = -dot(NM,CM)

y =linspace(CM(2,1)-10,CM(2,1)+10,51);
x =linspace(CM(1,1)-10,CM(1,1)+10,51);
[y2,x2]=meshgrid(y ,x );
z2=-(1/NM(3)).*(d+(NM(1).*x2)+(NM(2).*y2)); % ax+by+cz+d plane
figure (1)
% surf( x2,y2,z2 );
hold on

indpl1 = randi([0 length(x2)]);
indpl2 = randi([0 length(x2)]);

sym_plane = [CM,[x2(indpl1,indpl1),y2(indpl1,indpl1),z2(indpl1,indpl1)]',[x2(indpl2,indpl2),y2(indpl2,indpl2),z2(indpl2,indpl2)]']'




[um,un]=size(ivrs);

for fg=1:um
    %finding projection of upper points
    R = ivrs(fg,:)' - dot(ivrs(fg,:)' - CM, NM)*NM;
    
    p_ivrs(fg,:) =  2 * R - ivrs(fg,:)';
    
end


f=figure(1);
f.Position = [488 200 560 560];


%% plotting
figure(1)
% scatter3(vrs(:,1),vrs(:,2),vrs(:,3),5,'.b'); %right half
hold on

scatter3(CM(1,1),CM(2,1),CM(3,1),50,'*b'); %right half
scatter3(vrs(:,1),vrs(:,2),vrs(:,3),5,'.b'); %right half

hold on
% scatter3(ivrs(:,1),ivrs(:,2),ivrs(:,3),5,'or'); %right half

stlwrite(save_directory,fardad2.ConnectivityList, p_ivrs)

fardad3=stlread(save_directory);

scatter3(fardad3.Points(:,1),fardad3.Points(:,2),fardad3.Points(:,3),5,'og'); %right half
