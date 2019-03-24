%The goal here is to map one curve onto another by a rigid body
%transformation. Here, we optimization the transformation variables.

%Author: Abdullah Nazir

%Optimization variables: x, y, theta
clear all

load ('curves.mat')
load('energy_function_new.mat') %To plot quiver

x_bound = [0.5, 1];
y_bound = [0.0,0.5];
size = 15;

hx = (x_bound(2)-x_bound(1))/(size-1);
hy = (y_bound(2)-y_bound(1))/(size-1);
[U,V] = gradient(U_flex_matrix, hx, hy);

figure
axis equal
quiver(X,Y,U,V)
hold on


global oec_proc ftc_proc

ftc = transpose(finger_traj_curve);
oec = transpose(optimal_energy_curve);


[ftc_proc, oec_proc] = processCurves(ftc, oec);

x_init = 0.5*ones(1,3);

x_lb = [-2, -2, -pi];
x_ub = [2, 2, pi];


options = optimoptions('fmincon',...
                    'Display',...
                    'iter-detailed',...
                    'Algorithm',...
                    'interior-point');
               
x = fmincon(@objectiveFunction,...
            x_init,...
            [],[],[],[],x_lb,x_ub,...
            @constraintFunctions,...
            options)


ftc_trans = transformCurve(x(1),x(2),x(3), ftc_proc);


plot(ftc_trans(1,:), ftc_trans(2,:), 'LineWidth', 2.5)
hold on
plot(ftc_proc(1,:), ftc_proc(2,:), 'LineWidth', 2.5)
hold on
plot(oec_proc(1,:), oec_proc(2,:), 'LineWidth', 4, 'Color', 'g')


% lgd = legend('energy gradient field', ...
%        'transformed fingertip trajectory', ...
%        'untransformed fingertip trajectory', ...
%        'energy optimal path', 'location', 'northwest');
% 
% lgd.FontSize = 15;
% 
% xlabel('x')
% ylabel('y')
% title('Obtaining planar transform for fingertip trajectory to match energy optimal path')
% 



function fval = objectiveFunction(x)

    global oec_proc ftc_proc

    ftc_trans = transformCurve(x(1),x(2),x(3), ftc_proc);

    oec_proc_nonan = removeNaNs(oec_proc);
    ftc_proc_nonan = removeNaNs(ftc_trans);

    fval = distance_metric(ftc_proc_nonan, oec_proc_nonan);

end


function [c, c_eq] = constraintFunctions(x)

    c=[];
    c_eq=[];

end

function curve_nonan = removeNaNs(curve)

    curve_x = curve(1,:);
    curve_y = curve(2,:);

    curve_x_nonan = curve_x(~isnan(curve_y));
    curve_y_nonan = curve_y(~isnan(curve_y));

    curve_nonan = [curve_x_nonan; curve_y_nonan];
end


function [ftc_proc, oec_proc] = processCurves(ftc, oec)

    %compute bounds on x coordinate
    ul_bound_x = [min(ftc(1,:)), min(oec(1,:)), max(ftc(1,:)), max(oec(1,:))];
    
    x_new = linspace(min(ul_bound_x), max(ul_bound_x), 1000);
    
    
    ftc_y = interp1(ftc(1,:), ftc(2,:), x_new);
    ftc_proc = [x_new; ftc_y];
    
    
    oec_y = interp1(oec(1,:), oec(2,:), x_new);
    oec_proc = [x_new; oec_y];
        
end


function ftc_trans = transformCurve(x,y,theta, ftc)

    trans_mat = [cos(theta), -sin(theta), x;...
                 sin(theta), cos(theta), y;...
                 0,0,1];
    
    ftc_trans = zeros(2, size(ftc,2));
    
    for iter = 1:size(ftc,2)
        
            
        pt_trans = trans_mat*[ftc(1,iter); ftc(2,iter); 1];

        ftc_trans(:, iter) = [pt_trans(1); pt_trans(2)];
               
    end
    
end

function d_value = distance_metric(ftc_trans, oec_proc)

    distance_vec = zeros(1,size(ftc_trans,2)-size(oec_proc,2));
    
    oec_segment = oec_proc;
    
    
    for iter1 = 1:(size(ftc_trans,2)-size(oec_proc,2))
        
        ftc_segment = ftc_trans(:,iter1:iter1+size(oec_proc,2)-1);
        
        
        error_mat_squared = (oec_segment-ftc_segment).^2;
        
        
        distance_vec(iter1) = sum(sqrt(sum(error_mat_squared,1)));
             
           
           
    end
    
    d_value = min(distance_vec);

end

