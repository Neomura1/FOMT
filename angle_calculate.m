%% select interval
close all
clear
select_interval = 1; 
start_time = 322.5; end_time = 327.5; 
fit = 1; 
jump_time = 324.2; 
%% read data
[file_name,filefolder]=uigetfile({'*.tdms*'},'选择文件');
Data_Path=fullfile(filefolder,file_name);

tic
data = TDMS_getStruct(Data_Path);
temp_cell_ex =  fieldnames(data);

DNA_x = data.(temp_cell_ex{2}).Pos_X.data;
DNA_y = data.(temp_cell_ex{2}).Pos_Y.data;
DNA_z = data.(temp_cell_ex{2}).Pos_Z.data;
t = data.(temp_cell_ex{2}).Time.data;

N_cut = min([length(DNA_x),length(DNA_y),length(DNA_z),length(t)]);
DNA_x = data.(temp_cell_ex{2}).Pos_X.data(1:N_cut);
DNA_y = data.(temp_cell_ex{2}).Pos_Y.data(1:N_cut);
DNA_z = data.(temp_cell_ex{2}).Pos_Z.data(1:N_cut);
t = data.(temp_cell_ex{2}).Time.data(1:N_cut);
clear temp_cell_ex data

if select_interval == 1
    fs = length(t)/t(end);
    start_in = ceil(fs*start_time*60); end_in = floor(fs*end_time*60);
    if start_in <= 0
        start_in = 1;
    end
    if end_in > length(DNA_z)
        end_in = length(DNA_z);
    end
else
    start_in = 1; end_in = length(DNA_z);
end
DNA_x = DNA_x(start_in:end_in); DNA_y = DNA_y(start_in:end_in); DNA_z = DNA_z(start_in:end_in);
t = t(start_in:end_in)./60;

if fit == 1 & jump_time >= 0 & jump_time>start_time & jump_time<end_time
    jump_N = ceil((jump_time - start_time)*fs*60);
else
    jump_N = -1;
end
%% calculate angle
[xc,yc,R,~] = circfit(DNA_x,DNA_y);
clear beta;
data_number=length(DNA_x);
data_analysis_x=DNA_x-xc;
data_analysis_y=DNA_y-yc;

distance_2_center=sqrt((data_analysis_x).^2+(data_analysis_y).^2);
data_normal_x=data_analysis_x./distance_2_center;
data_normal_y=data_analysis_y./distance_2_center;

beta_count=1;
data_delete_x(1)=data_normal_x(1);
data_delete_y(1)=data_normal_y(1);
start_point=1;
end_point=2;
for i=1:data_number
    if end_point<=data_number
        distance_neighbour=sqrt((data_analysis_x(end_point)-data_analysis_x(start_point)).^2+(data_analysis_y(end_point)-data_analysis_y(start_point)).^2);
        if(distance_neighbour<0.98)
            beta_count=beta_count+1;
            data_delete_x(beta_count)=data_normal_x(end_point);
            data_delete_y(beta_count)=data_normal_y(end_point);
            start_point=end_point;
            end_point=start_point+1;
        else
            end_point=end_point+1;
        end
    end
end

data_delete_x=data_normal_x;
data_delete_y=data_normal_y;
data_number2=length(data_delete_x);
beta(1)=0;
for j=2:data_number2
    vector_0_x=data_delete_x(j);
    vector_0_y=data_delete_y(j);
    vector_1_x=data_delete_x(j-1);
    vector_1_y=data_delete_y(j-1);
    vector_0=[vector_0_x vector_0_y 0];
    vector_1=[vector_1_x vector_1_y 0];
    cross_vector=cross(vector_0,vector_1);
    beta_accelerate=sign(cross_vector(3))*acos(dot(vector_0,vector_1));
    if find(isnan(beta_accelerate)==1)
        beta_accelerate=0;
    end
    if abs(beta_accelerate)>10.62
        beta(j)=0;
    else
        beta(j)=beta(j-1)+real(beta_accelerate);
    end
end
turns = beta./2./pi;

%% plot
% xy plane
circle_plot(data_analysis_x,data_analysis_y,R,turns);
% R、turns、Ext.
cylindrical_plot(distance_2_center,turns,DNA_z,t,fit,jump_N);
disp(file_name);

toc
%% functions
function [xc,yc,R,a] = circfit(x,y)
n=length(x); xx=x.*x; yy=y.*y; xy=x.*y;
A=[sum(x) sum(y) n;sum(xy) sum(yy)...
sum(y);sum(xx) sum(xy) sum(x)];
B=[-sum(xx+yy) ; -sum(xx.*y+yy.*y) ; -sum(xx.*x+xy.*y)];
a=A\B;
xc = -.5*a(1);
yc = -.5*a(2);
R = sqrt((a(1)^2+a(2)^2)/4-a(3));
end
function circle_plot(x,y,r,turns)
% x&y
figure('name','circle');
s_temp = get(0);
w = s_temp.ScreenSize(3); h = s_temp.ScreenSize(4);
set(gcf,'unit','pixels','position',[0.01*w,0.6*h,0.6*h,0.3*h])

subplot(1,2,1);
scatter(x,y,3,'filled');
hold on
plot(0,0,'+');
hold on
alpha=0:0.1:2*pi;
circle_x=r*cos(alpha);
circle_y=r*sin(alpha);
plot(circle_x,circle_y,'g.');
axis equal;

xlabel('x(\mum)');ylabel('y(\mum)');
hold off
% angle
subplot(1,2,2);
turns_temp = turns - min(turns);
beta_normalize = 360*(turns_temp - floor(turns_temp)) -180;
histogram(beta_normalize,'Orientation','vertical','Normalization','pdf');
title('Angle/°');
end
function cylindrical_plot(distance_2_center,turns,DNA_z,time_calculate,fit,jump_N)
figure('name','R、\theta、z');
s_temp = get(0);
w = s_temp.ScreenSize(3); h = s_temp.ScreenSize(4);
set(gcf,'unit','pixels','position',[0.65*h,0.3*h,1.1*h,0.6*h])

for k = 1:18
    ax(k) = subplot(3,6,k);
end
subplot(3,6,[1 5]);
ax1 = plot(time_calculate,distance_2_center);
ylabel('R(\mum)'); yl1 = ylim;
subplot(3,6,[7 11]);
ax3 = plot(time_calculate,turns);
ylabel('turns'); yl3 = ylim;
subplot(3,6,[13 17]);
ax5 = plot(time_calculate,DNA_z);
xlabel('time(min)'),ylabel('Ext.(\mum)'); yl5 = ylim;

if fit == 1
    subplot(3,6,6);
    A1 = histogram(distance_2_center,'Orientation','horizontal','Normalization','probability','EdgeColor','auto'); hold on
    ylim(yl1); xlabel('Norm. Counts');
    options = statset('Display','off');
    gm1 = fitgmdist(distance_2_center',1,'Options',options);
    x = linspace(A1.BinLimits(1),A1.BinLimits(2),200);
    y = pdf(gm1,x');
    plot(y*A1.BinWidth,x,'Color','#A2142F','LineWidth',1.5)

    subplot(3,6,12);
    A2 = histogram(turns,'Orientation','horizontal','Normalization','probability','EdgeColor','auto'); hold on
    ylim(yl3); xlabel('Norm. Counts');
    gm2 = fitgmdist(turns',2,'Options',statset('Display','off','MaxIter',1500,'TolFun',1e-6));
    x = linspace(A2.BinLimits(1),A2.BinLimits(2),200);
    y1 = A2.BinWidth*gm2.ComponentProportion(1)*pdf('Normal',x, gm2.mu(1),sqrt(gm2.Sigma(1))); 
    y2 = A2.BinWidth*gm2.ComponentProportion(2)*pdf('Normal',x, gm2.mu(2),sqrt(gm2.Sigma(2)));
    plot(y1,x,'Color','#D95319','LineWidth',1.5); hold on
    plot(y2,x,'Color','#D95319','LineWidth',1.5); hold on
    plot(y1+y2,x,'Color','#A2142F','LineWidth',1.5); hold on
    text(1.1*max(y1+y2),gm2.mu(1),num2str(gm2.mu(1),'%.2f'),'VerticalAlignment','middle');
    text(1.1*max(y1+y2),gm2.mu(2),num2str(gm2.mu(2),'%.2f'),'VerticalAlignment','middle');
    disp([gm2.mu(1) gm2.mu(2)]);
    clear x y1 y2

    subplot(3,6,18);
    A3 = histogram(DNA_z,'Orientation','horizontal','Normalization','probability','EdgeColor','auto'); hold on
    ylim(yl5); xlabel('Norm. Counts');
    x = linspace(A3.BinLimits(1),A3.BinLimits(2),200);
    if jump_N < 0
        gm3 = fitgmdist(DNA_z',2,'Options',statset('Display','off','MaxIter',1500,'TolFun',1e-10));
        ComponentProportion = [gm3.ComponentProportion(1) gm3.ComponentProportion(2)];
        mu = [gm3.mu(1) gm3.mu(2)];
        sigma = [sqrt(gm3.Sigma(1)) sqrt(gm3.Sigma(2))];
    else
        gm3_1 = fitgmdist(DNA_z(1:jump_N)',1,'Options',statset('Display','off','MaxIter',1500,'TolFun',1e-10));
        gm3_2 = fitgmdist(DNA_z(jump_N:end)',1,'Options',statset('Display','off','MaxIter',1500,'TolFun',1e-10));
        N = length(DNA_z);
        ComponentProportion = [jump_N/N (N-jump_N)/N];
        mu = [gm3_1.mu gm3_2.mu];
        sigma = [sqrt(gm3_1.Sigma) sqrt(gm3_2.Sigma)];
    end
    y1 = A3.BinWidth*ComponentProportion(1)*pdf('Normal',x, mu(1), sigma(1));
    y2 = A3.BinWidth*ComponentProportion(2)*pdf('Normal',x, mu(2), sigma(2));
    plot(y1,x,'Color','#D95319','LineWidth',1.5); hold on
    plot(y2,x,'Color','#D95319','LineWidth',1.5); hold on
    plot(y1+y2,x,'Color','#A2142F','LineWidth',1.5); hold on
    text(1.1*max(y1+y2),mu(1),num2str(mu(1),'%.3f'),'VerticalAlignment','middle');
    text(1.1*max(y1+y2),mu(2),num2str(mu(2),'%.3f'),'VerticalAlignment','middle');
    disp([mu(1) mu(2)]);
else
    subplot(3,6,6);
    A1 = histogram(distance_2_center,'Orientation','horizontal','Normalization','probability','EdgeColor','auto'); hold on
    ylim(yl1); xlabel('Norm. Counts');

    subplot(3,6,12);
    A2 = histogram(turns,'Orientation','horizontal','Normalization','probability','EdgeColor','auto'); hold on
    ylim(yl3); xlabel('Norm. Counts');

    subplot(3,6,18);
    A3 = histogram(DNA_z,'Orientation','horizontal','Normalization','probability','EdgeColor','auto'); hold on
    ylim(yl5); xlabel('Norm. Counts');
end
end


