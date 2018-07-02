%% script for reading all measured signals

fprintf('\n ############################################################ \n')
fprintf(' ############ Starting to read all measured data ############ \n')
fprintf(' ############################################################ \n \n')

% folder ''Partikel Tracking 040816''
foldername = 'Partikel Tracking 040816';
fprintf(['\n reading folder **',foldername,'** ... \n'])
% all the subfolders listed in the 'dirOutput'
dirOutput=dir(fullfile(foldername,'*'));

%number of the folders
N = length(dirOutput);

% find the current path
path;
folder = cd;

%enter every single folder, and read all the tracks.txt
for i = 3:N
  fname = dirOutput(i).name;
  a = dir(fullfile(foldername,'/',fname,'*.txt'));
  cd(['./Partikel Tracking 040816/',fname]);
% cd('C:\Users\weizh\OneDrive\Desktop\Zoologie\Partikel Tracking 040816\040816 fish 4 film 6')
%   
  A = importdata('tracks.txt',' ');
  nLines = size(A);
for numLine = 2:nLines(1)
    B = cellstr(A{numLine,1});
    C = strread(B{1,1},'%s');
    nChar = size(C);
    for numChar = 1:nChar
        if ~isletter(C{numChar,1})
            MOut(numLine,numChar) = str2num(C{numChar,1});
        else
            MOut(numLine,numChar) = NaN;
        end
    end
end


ref = csvread('axis points.csv',1,5,[1,5,2,6]);

% 
% save(fname,'MOut');

%initiate all the numbers with tracking number
Velocity = MOut(:,1:6);
Velocity(:,6) = 0;

%find all the data with the same tracking number

Value = MOut(2,2);
NValue = find(MOut(:,2) == Value);
[XSize, YSize] = size(MOut);

%calculation of all the particle velocities
while NValue(end) < XSize
NValue = find(MOut(:,2)==Value);
Number_NValue = length(NValue);

for i_th_particle = NValue(1)+1 : NValue(end)

    Velocity(i_th_particle,6) =sqrt((MOut(i_th_particle,4) - MOut(i_th_particle-1,4))^2 +(MOut(i_th_particle,5) - MOut(i_th_particle-1,5))^2)/1.299;


    
end 




    %find the biggest velocity in all the velocities, showed in the 7th
    %column
    [Velocity_max, Pos_max] = max(Velocity(NValue(1):NValue(end), 6));
    Velocity(NValue(Pos_max),7) = Velocity_max; 
    
    %find the average velocity in all the velocities, showed in the 8th
    %column
    
    Velocity(NValue(end),8) = mean(Velocity(NValue(1):NValue(end), 6));
    
        
    %average velocity from the Start to the End, showed in the 9th column,
    %at the last of the tracking number
    Velocity(NValue(end),9) =sqrt((MOut(NValue(end),4) - MOut(NValue(1),4))^2 +(MOut(NValue(end),5) - MOut(NValue(1),5))^2)/(1.299*Number_NValue);

    
    %every single particle Runlengh: 10th column
    for i_th_particle = NValue(1)+1 : NValue(end)

    Velocity(i_th_particle,10) = Velocity(i_th_particle,6)*1.299;

    
    end
    %sum particle Runlength, 11th column
     Velocity(NValue(end),11) = sum(Velocity(NValue(1):NValue(end), 10));
    
    %find out wheather the particle goes in the direction a or b, every
    %single 
       %if the particle goes in the direction of ref_a  (oberer Punkt), in the
    %column 12, the last value of this tracking series is +1
    %otherwise, -1
    
for i_th_particle = NValue(1)+1 : NValue(end)

     %find out wheather the particle goes in the direction axis reference AB
    %Set the Reference 1, point A
    %Set the Reference 2, point B
    %The Particle ist temporarily C
    AB_x = ref(1,1) - ref(2,1);
    AB_y = ref(1,2) - ref(2,2);
    AB = sqrt(AB_x^2+AB_y^2);
    
    AC_x_old = ref(1,1) - MOut(i_th_particle-1,4);
    AC_y_old = ref(1,2) - MOut(i_th_particle-1,5);
    AC_old = sqrt((MOut(i_th_particle-1,4)  - ref(1,1))^2 +(MOut(i_th_particle-1,5) - ref(1,2))^2);
    
    AC_x_new = ref(1,1) - MOut(i_th_particle,4);
    AC_y_new = ref(1,2) - MOut(i_th_particle,5);
    AC_new = sqrt((MOut(i_th_particle,4)  - ref(1,1))^2 +(MOut(i_th_particle,5) - ref(1,2))^2);
    
    angle_old = acos((AB_x*AC_x_old + AB_y*AC_y_old)/AB/AC_old);
    distance_a_1 = abs(AC_old*sin(angle_old));
    
    angle_new = acos((AB_x*AC_x_new + AB_y*AC_y_new)/AB/AC_new);
    distance_b_1 = abs(AC_new*sin(angle_new));  
    
    if distance_a_1 < distance_b_1
        Velocity(i_th_particle,12) = 1;
    else
        Velocity(i_th_particle,12) = -1;
    end
end

     %find out wheather the particle goes in the direction axis reference AB
    AC_x_first = ref(1,1) - MOut(NValue(1),4);
    AC_y_first = ref(1,2) - MOut(NValue(1),5);
    AC_first = sqrt((MOut(NValue(1),4)  - ref(1,1))^2 +(MOut(NValue(1),5) - ref(1,2))^2);
    
    AC_x_last= ref(1,1) - MOut(i_th_particle,4);
    AC_y_last = ref(1,2) - MOut(i_th_particle,5);
    AC_last = sqrt((MOut(NValue(end),4)  - ref(1,1))^2 +(MOut(NValue(end),5) - ref(1,2))^2);
    
    angle_first = acos((AB_x*AC_x_first + AB_y*AC_y_first)/AB/AC_first);
    distance_a = abs(AC_first*sin(angle_first));
    
    angle_last = acos((AB_x*AC_x_last + AB_y*AC_y_last)/AB/AC_last);
    distance_b = abs(AC_last*sin(angle_last));  
    
    
    
    
    %if the particle goes in the direction of ref_a  (oberer Punkt), in the
    %column 13, the last value of this tracking series is +1
    %otherwise, -1
    if distance_a < distance_b
        Velocity(NValue(end),13) = 1;
    else
        Velocity(NValue(end),13) = -1;
    end 
    

     %sum of the values from the column 10 in column 14. Runlength.
    
    for i_th_particle = NValue(1) : NValue(end)

    Velocity(i_th_particle,14) =sum(Velocity(NValue(1):i_th_particle,10));


    end    
    
    
if NValue(end) < XSize
Value = MOut(NValue(end) + 1,2);

end 


filename = [fname,'_velocity.txt'];
%save all the velocities in the 6th column
save(filename, 'Velocity','-ascii');

     
end 
clear A B C a i i_th_particle N nChar nLines Number_NVlaue numChar numLine NValue Pos_max
clear MOut Value Velocity_max YSize

end


