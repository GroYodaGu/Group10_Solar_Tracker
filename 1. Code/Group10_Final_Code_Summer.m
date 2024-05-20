%% Group 5
% Mari Norwie
% Kashaka Sithole
% Christy Mkhari
% Ipeleng Tshabalala

% Code 100% Completed on 13/05/2024

syms s k;
stf = tf('s');

load sunPositionData.mat;
%% System Initialization
%%Shout-out to Matlab starting with Simulink that we can use their system.
% Motor parameters
Kf = 0.07;             % Back EMF cosntant, [V/(rad/s)]
Kt = 0.07;             % Torque constant, [N*m/A]
Kp = 0.318;              % Potentiometer gain,
L = 1e-5;              % Inductance, [H]
R = 10;                % Resistance, [Ohm]
Kg = 3000;                % Gear ratio, []
Td = 44.4;

% Solar panel parameters
m = 50;               % Mass, [kg]
w = 1.04;             % Width, [m]
l = 1.4;              % Length, [m]
d = 0.1;              % Depth, [m]
A = w*l;              % Area, [m^2]
Kd = 5;               % Damping constant, [N*m/(rad/s)]
Ap = 7;
sec = 10;
Lat = 26.2;
n = 355;
D = -23.45*sin((360/365)*(n+284));
local = 6.5;
h = (local-12)*15;
beta =90-(asin(sin(Lat)*sin(D)+cos(Lat)*cos(D)*cos(h))); % Elevation angle, [rad]
J = (m/12)*(l^2*cos(beta)^2+d^2*sin(beta)^2+w^2);
J0 = (m/12)*(l^2*cos(beta)^2+d^2*sin(beta)^2+w^2);

%% Numerators and Denominators Intialization

plant_num = Kg*Kt;
plant_den = (L*stf+R)*(J0*stf^2 + Kd*stf)+Kg^2*Kt*Kf*stf;

Motor_Num = 1;
Motor_Den = [L R];

Panel_Num = 1;
Panel_Den_lin = [J0 Kd];

%% Denominator Coefficients

A = 8.292e-05;
B = 82.92; 
C = 4.415e04; 
D = 21.24;
den_Open = [A B C 0];
den_Feedback = [A B C D];

%% Transfer Functions

plant_TF = plant_num/plant_den
Open_System_TF = Kp*plant_TF
Cl_system_TF = feedback(Kp*plant_TF,Kp)


%% Time-domain

figure(1); subplot(2,1,1); step(Open_System_TF); grid on; legend('Open-Loop System'); 
subplot(2,1,2); step(Cl_system_TF); grid on; legend('Feedback-Loop System'); 
stepinfo(Open_System_TF);
stepinfo(Cl_system_TF);

%% Frequency-domain

figure(2); bode(Open_System_TF); grid on; legend('Open-Loop System'); 
figure(3); bode(Cl_system_TF); grid on; legend('Feedback-Loop System');
figure(4); margin(Cl_system_TF);grid on;
[Gm,Pm,Wcg,Wcp] = margin(Cl_system_TF);

%% Frequency-domain Analytical Calculations 
r = B/A;
t = C/A;
roots(den_Feedback)

%% Routh=Hurwits

Routh_Open = Routh_Hurwitz(den_Open);
Routh_Feedback = Routh_Hurwitz(den_Feedback);

%% Nyquist Criterion

figure(4); nyquist(Open_System_TF); grid on; legend('Open-Loop System');
figure(5); nyquist(Cl_system_TF); grid on; legend('Feedback-Loop System');

%% Nyquist Analytical Analysis

W = -500:0.01:500;
Re = (D^2-D*B*W.^2)./(-A^2*W.^6+(2*C*A+B^2)*W.^4-(2*D*B+C^2)*W.^2+D^2);
Im = -(D*(C*W-A*W.^3))./(-A^2*W.^6+(2*C*A+B^2)*W.^4-(2*D*B+C^2)*W.^2+D^2);
% Magnitude = 1./sqrt((D-B*W.^2)^2-(C*W-A*W.^3)^2);
% Phase = -atan(D*(C*W-A*W.^3))*57.3;
figure;plot(Re,Im);
RE = (-Kg*Kt*B*W.^2)./(-A^2*W.^6+(2*C*A+B^2)*W.^4-C^2*W.^2);
IM = -(Kg*Kt*(C*W-A*W.^3))./(-A^2*W.^6+(2*C*A+B^2)*W.^4-C^2*W.^2);
figure;plot(RE,IM);

%% BIBO Stability

pole_Open = pole(Open_System_TF);
pole_Feedback = pole(Cl_system_TF);

%% State-Space

Linear_system = ss(A,B,C,D);
pole_sys = pole(Linear_system);


%% PID Selection
%% Proportional Control
t = 0:0.5:10;
KP = 100;
C_P = pid(KP);
sys_feedback_P1 = feedback(Kp*C_P*plant_TF,Kp);
step(sys_feedback_P1);
stepinfo(sys_feedback_P1);

KP = 500;
C_P = pid(KP);
sys_feedback_P2 = feedback(Kp*C_P*plant_TF,Kp);
step(sys_feedback_P2);
stepinfo(sys_feedback_P2);

KP = 1000;
C_P = pid(KP);
sys_feedback_P3 = feedback(Kp*C_P*plant_TF,Kp);
step(sys_feedback_P3);
stepinfo(sys_feedback_P3);

KP = 2000;
C_P = pid(KP);
sys_feedback_P4 = feedback(Kp*C_P*plant_TF,Kp);
step(sys_feedback_P4);
stepinfo(sys_feedback_P4);

KP = 3000;
C_P = pid(KP);
sys_feedback_P5 = feedback(Kp*C_P*plant_TF,Kp);
step(sys_feedback_P5);
stepinfo(sys_feedback_P5);

KP = 4000;
C_P = pid(KP);
sys_feedback_P6 = feedback(Kp*C_P*plant_TF,Kp);
step(sys_feedback_P6);
stepinfo(sys_feedback_P6);

KP = 4500;
C_P = pid(KP);
sys_feedback_P7 = feedback(Kp*C_P*plant_TF,Kp);
step(sys_feedback_P7);
stepinfo(sys_feedback_P7);

KP = 4700;
C_P = pid(KP);
sys_feedback_P8 = feedback(Kp*C_P*plant_TF,Kp);
step(sys_feedback_P8);
stepinfo(sys_feedback_P8);

figure(6); subplot(2,2,1); step(sys_feedback_P1,t);hold on; stepplot(sys_feedback_P2,t,'r--'); grid on; legend('Kp = 100', 'Kp = 500');
subplot (2,2,2);step(sys_feedback_P3,t); hold on; stepplot(sys_feedback_P4,t,'r--'); grid on; legend('Kp = 1000','Kp = 2000');
subplot (2,2,3);step(sys_feedback_P5,t); hold on; stepplot(sys_feedback_P6,t,'r--'); grid on; legend('Kp = 3000','Kp = 4000');
subplot (2,2,4);step(sys_feedback_P7,t); hold on; stepplot(sys_feedback_P8,t,'r--'); grid on; legend('Kp = 4500','Kp = 4700');

%% Proportional & Integral Control
KP = 240;
KI = 180;
C_PI = pid(KP,KI);
sys_feedback_PI1 = feedback(Kp*C_PI*plant_TF, Kp);
step(sys_feedback_PI1);
stepinfo(sys_feedback_PI1);

KP = 500;
KI = 180;
C_PI = pid(KP,KI);
sys_feedback_PI2 = feedback(Kp*C_PI*plant_TF, Kp);
step(sys_feedback_PI2);
stepinfo(sys_feedback_PI2);

KP = 500;
KI = 50;
C_PI = pid(KP,KI);
sys_feedback_PI3 = feedback(Kp*C_PI*plant_TF, Kp);
step(sys_feedback_PI3);
stepinfo(sys_feedback_PI3);

KP = 900;
KI = 50;
C_PI = pid(KP,KI);
sys_feedback_PI4 = feedback(Kp*C_PI*plant_TF, Kp);
step(sys_feedback_PI4);
stepinfo(sys_feedback_PI4);

KP = 1000;
KI = 50;
C_PI = pid(KP,KI);
sys_feedback_PI5 = feedback(Kp*C_PI*plant_TF, Kp);
step(sys_feedback_PI5);
stepinfo(sys_feedback_PI5);

KP = 1500;
KI = 50;
C_PI = pid(KP,KI);
sys_feedback_PI6 = feedback(Kp*C_PI*plant_TF, Kp);
step(sys_feedback_PI6);
stepinfo(sys_feedback_PI6);

KP = 3000;
KI = 100;
C_PI = pid(KP,KI);
sys_feedback_PI7 = feedback(Kp*C_PI*plant_TF, Kp);
step(sys_feedback_PI7);
stepinfo(sys_feedback_PI7);

KP = 4700;
KI = 100;
C_PI = pid(KP,KI);
sys_feedback_PI8 = feedback(Kp*C_PI*plant_TF, Kp);
step(sys_feedback_PI8);
stepinfo(sys_feedback_PI8);


figure(7); subplot(2,2,1); step(sys_feedback_PI1,t);hold on; stepplot(sys_feedback_PI2,t,'r--'); grid on; legend('Kp = 240 Ki = 180', 'Kp = 500 Ki = 180');
subplot (2,2,2); step(sys_feedback_PI3,t); hold on; stepplot(sys_feedback_PI4,t,'r--'); grid on; legend('Kp = 500 Ki = 50','Kp = 900 Ki = 50');
subplot (2,2,3); step(sys_feedback_PI5,t); hold on; stepplot(sys_feedback_PI6,t,'r--'); grid on; legend('Kp = 1000 Ki = 50','Kp = 1500 Ki = 50');
subplot (2,2,4); step(sys_feedback_PI7,t); hold on; stepplot(sys_feedback_PI8,t,'r--'); grid on; legend('Kp = 3000 Ki = 100','Kp = 4700 Ki = 100');

%% Proportional & Derivative Control
KP = 500;
KI = 0;
KD = 100;
C_PD = pid(KP,KI,KD);
sys_feedback_PD1 = feedback(Kp*C_PD*plant_TF,Kp);
step(sys_feedback_PD1);
stepinfo(sys_feedback_PD1);

KP = 700;
KI = 0;
KD = 100;
C_PD = pid(KP,KI,KD);
sys_feedback_PD2 = feedback(Kp*C_PD*plant_TF,Kp);
step(sys_feedback_PD2);
stepinfo(sys_feedback_PD2);

KP = 900;
KI = 0;
KD = 100;
C_PD = pid(KP,KI,KD);
sys_feedback_PD3 = feedback(Kp*C_PD*plant_TF,Kp);
step(sys_feedback_PD3);
stepinfo(sys_feedback_PD3);

KP = 4700;
KI = 0;
KD = 100;
C_PD = pid(KP,KI,KD);
sys_feedback_PD4 = feedback(Kp*C_PD*plant_TF,Kp);
step(sys_feedback_PD4);
stepinfo(sys_feedback_PD4);

KP = 4700;
KI = 0;
KD = 50;
C_PD = pid(KP,KI,KD);
sys_feedback_PD5 = feedback(Kp*C_PD*plant_TF,Kp);
step(sys_feedback_PD5);
stepinfo(sys_feedback_PD5);

KP = 4700;
KI = 0;
KD = 20;
C_PD = pid(KP,KI,KD);
sys_feedback_PD6 = feedback(Kp*C_PD*plant_TF,Kp);
step(sys_feedback_PD6);
stepinfo(sys_feedback_PD6);

KP = 4700;
KI = 0;
KD = 10;
C_PD = pid(KP,KI,KD);
sys_feedback_PD7 = feedback(Kp*C_PD*plant_TF,Kp);
step(sys_feedback_PD7);
stepinfo(sys_feedback_PD7);

KP = 4700;
KI = 0;
KD = 5;
C_PD = pid(KP,KI,KD);
sys_feedback_PD8 = feedback(Kp*C_PD*plant_TF,Kp);
step(sys_feedback_PD8);
stepinfo(sys_feedback_PD8);

figure(8); subplot(2,2,1); step(sys_feedback_PD1,t);hold on; stepplot(sys_feedback_PD2,t,'r--'); grid on; legend('Kp = 500 Kd = 100', 'Kp = 700 Kd = 100');
subplot (2,2,2); step(sys_feedback_PD3,t); hold on; stepplot(sys_feedback_PD4,t,'r--'); grid on; legend('Kp = 900 Kd = 100','Kp = 4700 Kd = 100');
subplot (2,2,3); step(sys_feedback_PD5,t); hold on; stepplot(sys_feedback_PD6,t,'r--'); grid on; legend('Kp = 4700 Kd = 100','Kp = 4700 Kd = 20');
subplot (2,2,4); step(sys_feedback_PD7,t); hold on; stepplot(sys_feedback_PD8,t,'r--'); grid on; legend('Kp = 4700 Kd = 10','Kp = 4700 Kd = 5');

%% Proportional, Integral and Derivative Control
KP = 500;
KI = 40;
KD = 100;
C_PID = pid(KP,KI,KD);
sys_feedback_PID1 = feedback(Kp*C_PID*plant_TF,Kp);
step(sys_feedback_PID1);
stepinfo(sys_feedback_PID1);

KP = 900;
KI = 40;
KD = 90;
C_PID = pid(KP,KI,KD);
sys_feedback_PID2 = feedback(Kp*C_PID*plant_TF,Kp);
step(sys_feedback_PID2);
stepinfo(sys_feedback_PID2);

KP = 1000;
KI = 50;
KD = 90;
C_PID = pid(KP,KI,KD);
sys_feedback_PID3 = feedback(Kp*C_PID*plant_TF,Kp);
step(sys_feedback_PID3);
stepinfo(sys_feedback_PID3);

KP = 3000;
KI = 50;
KD = 50;
C_PID = pid(KP,KI,KD);
sys_feedback_PID4 = feedback(Kp*C_PID*plant_TF,Kp);
step(sys_feedback_PID4);
stepinfo(sys_feedback_PID4);

KP = 3500;
KI = 70;
KD = 30;
C_PID = pid(KP,KI,KD);
sys_feedback_PID5 = feedback(Kp*C_PID*plant_TF,Kp);
step(sys_feedback_PID5);
stepinfo(sys_feedback_PID5);

KP = 4000;
KI = 80;
KD = 20;
C_PID = pid(KP,KI,KD);
sys_feedback_PID6 = feedback(Kp*C_PID*plant_TF,Kp);
step(sys_feedback_PID6);
stepinfo(sys_feedback_PID6);

KP = 4500;
KI = 90;
KD = 10;
C_PID = pid(KP,KI,KD);
sys_feedback_PID7 = feedback(Kp*C_PID*plant_TF,Kp);
step(sys_feedback_PID7);
stepinfo(sys_feedback_PID7);

KP = 4700;
KI = 100;
KD = 5;
C_PID = pid(KP,KI,KD);
sys_feedback_PID8 = feedback(Kp*C_PID*plant_TF,Kp);
step(sys_feedback_PID8);
stepinfo(sys_feedback_PID8);

figure(9); subplot(2,2,1); step(sys_feedback_PID1,t);hold on; stepplot(sys_feedback_PID2,t,'r--'); grid on; legend('Kp = 500 Ki = 40 Kd = 100', 'Kp = 700 Ki = 40 Kd = 90');
subplot (2,2,2);step(sys_feedback_PID3,t); hold on; stepplot(sys_feedback_PID4,t,'r--'); grid on; legend('Kp = 1000 Ki = 50 Kd = 90','Kp = 3000 Ki = 50 Kd = 50');
subplot (2,2,3);step(sys_feedback_PID5,t); hold on; stepplot(sys_feedback_PID6,t,'r--'); grid on; legend('Kp = 3500 Ki = 70 Kd = 30','Kp = 4000 Ki = 80 Kd = 20');
subplot (2,2,4);step(sys_feedback_PID7,t); hold on; stepplot(sys_feedback_PID8,t,'r--'); grid on; legend('Kp = 4500 Ki = 90 Kd = 10','Kp = 4700 Ki = 100 Kd = 5');

figure(10); subplot(2,2,1); step(sys_feedback_P8,t);grid on; legend('Kp = 4700');
subplot (2,2,2); step(sys_feedback_PI8,t); grid on; legend('Kp = 4700 Ki = 100');
subplot (2,2,3); step(sys_feedback_PD5,t); grid on; legend('Kp = 4700 Kd = 50');
subplot (2,2,4); step(sys_feedback_PID7,t); grid on; legend('Kp = 4500 Ki = 90 Kd = 10');
%% Root Locus

figure(11); rlocus(Cl_system_TF); legend('Feedback-Loop System');

%% Obtaining Simulink data

% Time Variable
time = out.tout;

% Linear Vs Non-linear with Various inputs
% StepLin = out.Linear_Step_Out;
% StepNLin = out.Non_Linear_Step_Out;
% RampLin = out.Linear_Ramp_Out;
% RampNLin = out.Non_Linear_Ramp_Out;
% SineLin = out.Linear_Sine_Out;
% SineNLin = out.Non_Linear_Sine_Out;
% figure(12); plot(time,StepLin,'g'); hold on; plot(time,StepNLin,'r--'); grid on; title('Step Input Linear vs Non-Linear'); legend('Linear','Non-Linear'); xlabel('Time (s)'); ylabel('Amplitude');
% figure(13); plot(time,RampLin,'g'); hold on; plot(time,RampNLin,'b--'); grid on; title('Ramp Input Linear vs Non-Linear'); legend('Linear','Non-Linear'); xlabel('Time (s)'); ylabel('Amplitude');
% figure(14); plot(time,SineLin,'g'); hold on; plot(time,SineNLin,'m--'); grid on; title('Sine Input Linear vs Non-Linear'); legend('Linear','Non-Linear'); xlabel('Time (s)'); ylabel('Amplitude');

% Absolute Stability
% Pulse = out.Abs_Stability_feedback_pulse;
% figure(15); plot(time, Pulse); title('Absolute Stability with Pulse Wave'); grid on;

% BIBO Stability
% SineW = out.BIBO_Stability_feedback_Sin;
% Square = out.BIBO_Stability_feedback_Square;
% figure(16); plot(time,SineW); grid on; title('BIBO Stability with Sine Wave Input')  
% figure(17); plot(time,Square); grid on; title('BIBO Stability with Square Wave Input');

% System Behaviour with P Controller - No Disturbance
% AngleIn = out.Linear_PID_In_B;
% AngleOut = out.Linear_PID_OUT_B;
% figure(18); plot(time, AngleIn,'g'); hold on; plot(time, AngleOut,'r--'); grid on; title('Feedback System Behaviour with P Contoller - No Disturbance - Summer Time');
% legend('Sun Position','Panel Tracking'); xlabel('Time (s)'); ylabel('Sun Position (rad)');

% System Behaviour with P Controller - No Disturbance - Non-Linear
% AngleIn = out.Non_Linear_PID_In_B;
% AngleOut = out.Non_Linear_PID_OUT_B;
% figure(19); plot(time, AngleIn,'g'); hold on; plot(time, AngleOut,'r--'); grid on; title('Non-Linear Feedback System Behaviour with P Contoller - No Disturbance - Summer Time');
% legend('Sun Position','Panel Tracking'); xlabel('Time (s)'); ylabel('Sun Position (rad)');

% System Behaviour with Constant Disturbance
% AngleIn = out.Linear_PID_In_ConstD;
% AngleOut = out.Linear_PID_Out_ConstD;
% figure(20); plot(time,AngleIn,'g'); hold on; plot(time, AngleOut,'m--'); grid on; title('Feedback System Behaviour with P Controller and Constant Disturbance - Summer Time');
% legend('Sun Position','Panel Tracking'); xlabel('Time (s)'); ylabel('Sun Position (rad)');

% System Behaviour with Pulse Disturbance
% AngleIn = out.Linear_PID_In_PulseD;
% AngleOut = out.Linear_PID_Out_PulseD;
% figure(21); plot(time,AngleIn,'g'); hold on; plot(time, AngleOut,'b--'); grid on; title('Feedback System Behaviour with P Controller and Pulse Disturbance - Summer Time');
% legend('Sun Position','Panel Tracking'); xlabel('Time (s)'); ylabel('Sun Position (rad)');

% System Behaviour with Sine Wave Disturbance
% AngleIn = out.Linear_PID_In_SineD;
% AngleOut = out.Linear_PID_Out_SineD;
% figure(22); plot(time,AngleIn,'g'); hold on; plot(time, AngleOut,'k--'); grid on; title('Feedback System Behaviour with P Controller and Sine Wave Disturbance - Summer Time');
% legend('Sun Position','Panel Tracking'); xlabel('Time (s)'); ylabel('Sun Position (rad)');

%% Routh-Hurwitz Table Function
%%Shoutout to Joe Schenkel for this Routh-Hurwitz Table code that works perfectly
function [routh_table] = Routh_Hurwitz(Den)

syms e;             %create a symbol for epsilon (in case zero's show up in Routh Table)
syms k1;            %create symbol for K (in case user enters that value)
syms k2;
eps_flag = 0;       %flag for if a first column zero appears (e replaces it)
row_of_zeros = 0;   %variable for the power of s in which a row of zeros occured
constant_flag = 0;  %flag for if user is using symbolic constants

den_poly = Den;
fprintf('\r\n');

%%%%%%%%%%%%%%%%%%% Create the first 2 rows from input information %%%%%%%%%%%%%%
column_size = length(den_poly);
routh_table_1 = den_poly(1);
for x = 3 : 2 : column_size
    routh_table_1 = [routh_table_1, den_poly(x)];
end

routh_table_2 = den_poly(2);


for x = 4 : 2 : column_size
    routh_table_2 = [routh_table_2, den_poly(x)];
end

if(rem(column_size, 2) == 1)  % even polynomial
    routh_table_2 = [routh_table_2, 0];
    new_column_size = ((column_size + 1) / 2) - 1;  % 3rd column size for even polynomials
else
    new_column_size = (column_size / 2) - 1;  % 3rd column size for odd polynomials  
end

routh_table = [routh_table_1; routh_table_2];


%%%%%%%%%%%%%%%%%%% First 2 rows created, now finish the table %%%%%%%%%%%%%%

for z = 1 : 1 : column_size - 2  % each loop is a new row
    determinant = [routh_table(z,1) routh_table(z,2); routh_table(z+1, 1) routh_table(z+1,2)]; 
    routh_table_temp = - (det(determinant) / routh_table(z+1,1));

    if routh_table_temp == 0    %% Test to see if zero will be in table
        routh_table_temp = e;   %% Replace the zero with epsilon
        eps_flag = 1;
    end
    routh_table_new_row = routh_table_temp;

    for x = 2 : 1 : new_column_size     % Each loop is for a column in the table
        determinant = [routh_table(z,1) routh_table(z,x+1); routh_table(z+1,1) routh_table(z+1,x+1)]; 
        routh_table_temp = -(det(determinant) / routh_table(z+1,1));
        routh_table_new_row = [routh_table_new_row, routh_table_temp];
    end
    
    for y = length(routh_table_new_row) : 1 : new_column_size  % Fill the last columns with zeros
        routh_table_new_row = [routh_table_new_row, 0];
    end
    
    if row_of_zeros == 0   %% These if statements are for determining if a row of zeros has appeared
        if eps_flag == 1
            if routh_table_new_row(1, 2:new_column_size) == zeros(1,(new_column_size-1))
                row_of_zeros = column_size - 2 - z;
            end
        end
    end
    routh_table = [routh_table; routh_table_new_row];  % Add the new row to the table
end

%%%%%%%%%%%% Now it's time to display the Routh Table (and warnings) %%%%%%%%%%%
fprintf('\r\n');
fprintf(' Routh Table:\r\n\r\n');
if row_of_zeros > 0  % a row of zeros occurred, alert user of where and what to do
    disp(' This means you might have poles on the jw-axis');    
    fprintf('\r\n');
    disp(' To complete the table, please do the following (if necessary):');
    disp(sprintf('\t1.  Please differentiate the s^%d row', (row_of_zeros+1)));
    disp(sprintf('\t2.  Substitute answer in s^%d row', row_of_zeros));
    disp(sprintf('\t3.  Re-run program with s^%d as the highest order of the polynomial', (row_of_zeros+1)));
    disp(sprintf('\t    and the differentiated answer as the s^%d row.  \r\n\r\n', (row_of_zeros)));    
end

sym(routh_table);
disp(sprintf(' Highest order of Polynomial Entered = s^%d',(column_size-1)));

if constant_flag == 1       %Use pretty print if constants used
    pretty(routh_table);    
elseif eps_flag == 1        %Use pretty print if symbol 'e' used
    pretty(routh_table);
else                        %Use regular print otherwise
    disp(routh_table);
end
fprintf('\r\n\r\n\r\n');

end