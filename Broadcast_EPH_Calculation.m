% Jason Lee <jason_lee@sjtu.edu.cn>
% SJTU SEIEE
% Created: November 8, 2023

clear all
clc
close all

% Load Data
data_read = fopen('.\data_process.rnx', 'r');

if data_read == -1
    error('Failed openning file %s', 'D:\projects\renix4_error_calc\data_process.rnx');
end

pattern = '\d{4}|\d\d|-?\d\.\d{12}e[+-]\d\d';
parameters = [];

while true
    % Read file line by line
    line = fgetl(data_read);
    % Check if ends up the file or not
    if line == -1
        break;
    end
    
    % Extract values needed
    data = regexp(line, pattern, 'match');
    
    % Check if there are value needed
    if ~isempty(data)
        % Add the value to parameter array
        parameters{end+1, 1} = data;
    end
end


% Close file
fclose(data_read);

% Constant variables
mu = 3.986005e14 % GM
omega_e = 7.292115147e-5
PI = 3.1415926535898
c = 299792458

% Read Data
t_oe = str2double(parameters{5}{1})   %5,1
t_tm = str2double(parameters{11}{1})   %11,1
sqrt_a = str2double(parameters{4}{4}) %4,4
DELTA_n = str2double(parameters{3}{3})  %3,3
M_0 = str2double(parameters{3}{4})    %3,4
e = str2double(parameters{4}{2}) %4,2
a_bias = str2double(parameters{2}{8}) %2,4(2,8 after fixing)
a_drift = str2double(parameters{2}{9}) %2,5(2,9 after fixing)
a_dr = str2double(parameters{2}{10}) %2,6(2,10 after fixing)
i_0 = str2double(parameters{6}{1}) %6,1
OMEGA_0 = str2double(parameters{5}{3}) %5,3
omega = str2double(parameters{6}{3}) %6,3
C_rc = str2double(parameters{6}{2}) %6,2
C_uc = str2double(parameters{4}{1}) %4,1
C_us = str2double(parameters{4}{3}) %4,3
C_rs = str2double(parameters{3}{2}) %3,2
C_ic = str2double(parameters{5}{2}) %5,2
C_is = str2double(parameters{5}{4}) %5,4
OMEGA_dot = str2double(parameters{6}{4}) %6,4; Diffrential of OMEGA
i_dot = str2double(parameters{7}{1}) %7,1; Diffrential of i

a = sqrt_a^2;

% Step 1: Compute mean anomaly and mean motion
% Compute the relative time, where t is the time of transmission
% t_k = t - t_oe
t_k = 30
% Compute the mean motion
n_o = sqrt(mu/(a^3))
% Correct the mean motion
n = n_o + DELTA_n
% Compute the anomaly at time of transmission
M_k = M_0 + n*t_k

% Step 2: Solve eccentric anomaly
% Iteratively solve for eccentric anomaly
% 3 times iteration
E_k = M_k;
for i = 1:100
    E_k = M_k + e*sin(E_k);
end

% Step 3: Compute true anomaly
v_k = atan((sqrt(1-e^2)*sin(E_k))/(cos(E_k)-e))

% Step 4: Compute argument of latitude (of the satellite)
PHI_k = omega + v_k

% Step 5: Compute corrections to Keplerian orbit
delta_u_k = C_uc*cos(2*PHI_k) + C_us*sin(2*PHI_k)
delta_r_k = C_rc*cos(2*PHI_k) + C_rs*sin(2*PHI_k)
delta_i_k = C_ic*cos(2*PHI_k) + C_is*sin(2*PHI_k)

% Step 6: Compute corrected values
u_k = PHI_k + delta_u_k
r_k = a*(1 - e*cos(E_k)) + delta_r_k
i_k = i_0 + i_dot*t_k + delta_i_k

% Step 7: Correct longitude of ascending node
OMEGA_k = (OMEGA_0 - omega_e*t_oe) + (OMEGA_dot - omega_e)*t_k

% Step 8: Compute position at time of transmission
x = r_k*cos(u_k)
y = r_k*sin(u_k)

X_k = x*cos(OMEGA_k) - y*cos(i_k)*sin(OMEGA_k)
Y_k = x*sin(OMEGA_k) + y*cos(i_k)*cos(OMEGA_k)
Z_k = y*sin(i_k)

% Clock bias calculation
FF=-2*sqrt(mu)/(c^2);
DetTr = FF*e*sqrt_a*sin(E_k);
Clock_bias = a_bias +a_drift*t_oe+a_dr*(t_oe^2)+DetTr;

% result = [
%             cos(-OMEGA_k),-sin(-OMEGA_k),0;
%             sin(-OMEGA_k),cos(-OMEGA_k),0;
%             0,0,1
%         ] * [
%             1,0,0;
%             0,cos(-i_k),-sin(-i_k);
%             0,sin(-i_k),cos(-i_k)
%         ] * [
%             cos(-u_k),-sin(-u_k),0;
%             sin(-u_k),cos(-u_k),0;
%             0,0,1
%         ] * ([
%             r_k,0,0
%         ].')

% Step 9: Account for Earth rotation during propagation time
% result_p = [
%             cos(omega_e*t_tm),-sin(omega_e*t_tm),0;
%             sin(omega_e*t_tm),cos(omega_e*t_tm),0;
%             0,0,1
%         ]*result

% Save the result values
% Defines the variable to be stored and the corresponding variable name
variables = {
    '长半轴a', a;
'平均角速度n0', n_o;
'观测历元到参考历元的时间差tk', t_k;
'改正平均角速度n', n;
'平近点角Mk', M_k;
'偏近点角Ek', E_k;
'真近点角vk', v_k;
'纬度幅角Φk', PHI_k;
'纬度幅角改正项δuk', delta_u_k;
'径向改正项δrk', delta_r_k;
'轨道倾角改正项δik', delta_i_k;
'改正后的纬度幅角uk', u_k;
'改正后的径向rk', r_k;
'改正后的轨道倾角ik', i_k;
'卫星在轨道平面内的坐标x', x;
'卫星在轨道平面内的坐标y', y;
'历元升交点经度Xk', X_k;
'历元升交点经度Yk', Y_k;
'历元升交点经度Zk', Z_k;
'卫星钟差', Clock_bias;
};

filename = 'variables.txt';

% Open the file for writing
fileID = fopen(filename, 'w');

% Loop writes variable names and values to the file
for i = 1:size(variables, 1)
    varName = variables{i, 1};
    varValue = variables{i, 2};
    fprintf(fileID, '%s: %e\n', varName, varValue);
end

% Close file
fclose(fileID);
disp('Values are already saved in variable.txt under the main code path')

