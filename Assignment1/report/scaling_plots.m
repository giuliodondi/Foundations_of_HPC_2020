clc;close all;
format long;
figure_size=1024;


cores=[1,4,8,12,16,20,24,28,32,36,40,44,48];
runs=[1e+8;1e+9;1e+10;1e+11];
N_vec = logspace(8,11,100);

serial_time_8(:)=strong_ext08(1,:);
serial_time_9(:)=strong_ext09(1,:);
serial_time_10(:)=strong_ext10(1,:);
serial_time_11(:)=strong_ext11(1,:);

svs1p = zeros(4,2);

svs1p(1,:) = scaling_err(strong_ext08(2,:),serial_time_8);
svs1p(2,:) = scaling_err(strong_ext09(2,:),serial_time_9);
svs1p(3,:) = scaling_err(strong_ext10(2,:),serial_time_10);
svs1p(4,:) = scaling_err(strong_ext11(2,:),serial_time_11);

mpi_init_t = 1 + (strong_ext08(2,1) - serial_time_8(1,1))*(1e+08)/( serial_time_8(1,1) )*1./N_vec ;

figure(1)
hold all 

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);

title('Serial vs 1-core Parallel')
xlabel('Problem Size (N)');
ylabel('T_p(1) / T_s');
xlim([0.8e+8,1.2e+11]);
ylim([0.998,1.125]);
%xticks([1,2,3,4])
%xticklabels({'10^8','10^9','10^{10}','10^{11}'})

oneoverN_f = @(a,N) 1 + a*1e08./N;
oneoverN0 = [0.3];

oneoverN = lsqcurvefit(oneoverN_f, oneoverN0, runs ,svs1p(:,1) )

(strong_ext08(2,1) - serial_time_8(1,1))/( serial_time_8(1,1) )

errorbar(runs,svs1p(:,1),svs1p(:,2),'*','LineStyle','none','Color','k','LineWidth',1.5,'MarkerSize',3);
plot(N_vec,oneoverN_f(oneoverN,N_vec),'r','LineStyle',':','LineWidth',1.5);

set(gca,'xscale','log')


figure(11)
hold all 

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);

title('Serial time vs Problem Size')
xlabel('Problem Size (N)');
ylabel('T_s (seconds)');
xlim([0.8e+8,1.2e+11]);
%xticks([1,2,3,4])
%xticklabels({'10^8','10^9','10^{10}','10^{11}'})
serial_times=[serial_time_8(:),serial_time_9(:),serial_time_10(:),serial_time_11(:)];

line_f = @(c,x) transpose(c(1) + x.*c(2));
serial_coef0 = [0,1];

serial_coef = lsqcurvefit(line_f, serial_coef0, runs,serial_times(1,:),[0,0],[inf,inf])

serial_fit = line_f(serial_coef,runs); 

errorbar(runs,serial_times(1,:),serial_times(2,:),'*','LineStyle','none','Color','k','LineWidth',1.5,'MarkerSize',3);
plot(runs,serial_fit,'r','LineStyle',':','LineWidth',1.5);

set(gca,'xscale','log')
set(gca,'yscale','log')



ext_scaling_strong08 = scaling_err(serial_time_8,strong_ext08(2:end,:));
ext_scaling_strong09 = scaling_err(serial_time_9,strong_ext09(2:end,:));
ext_scaling_strong10 = scaling_err(serial_time_10,strong_ext10(2:end,:));
ext_scaling_strong11 = scaling_err(serial_time_11,strong_ext11(2:end,:));

int_scaling_strong08 = scaling_err(serial_time_8,strong_int08(:,:));
int_scaling_strong09 = scaling_err(serial_time_9,strong_int09(:,:,:));
int_scaling_strong10 = scaling_err(serial_time_10,strong_int10(:,:,:));
int_scaling_strong11 = scaling_err(serial_time_11,strong_int11(:,:,:));





figure(22)

hold all 
grid on
title('Strong Scaling, extrapolation, no overhead')
xlabel('Cores');
ylabel('Speedup S');

coress=0:1:100000;

 
th_scaling_noov(100000,serial_coef,1e08)

plot(coress,th_scaling_noov(coress,serial_coef,1e08) );
plot(coress,th_scaling_noov(coress,serial_coef,1e09) );
plot(coress,th_scaling_noov(coress,serial_coef,1e010) );
plot(coress,th_scaling_noov(coress,serial_coef,1e011) );






figure(2)

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
subplot(2,2,1)
hold all 
grid on
title('Strong Run-times 10-08')
xlabel('Cores');
ylabel('T_p (seconds');
xlim([1,48]);
xticks(cores);


errorbar(cores,strong_int08(:,1),strong_int08(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,strong_ext08(2:end,1),strong_ext08(2:end,2),'LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');


subplot(2,2,2)
hold all 
grid on
title('Strong Run-times 10-09')
xlabel('Cores');
ylabel('T_p (seconds');
xlim([1,48]);
xticks(cores);


errorbar(cores,strong_int09(:,1),strong_int09(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,strong_ext09(2:end,1),strong_ext09(2:end,2),'LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,3)
hold all 
hold all 
grid on
title('Strong Run-times 10-10')
xlabel('Cores');
ylabel('T_p (seconds');
xlim([1,48]);
xticks(cores);


errorbar(cores,strong_int10(:,1),strong_int10(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,strong_ext10(2:end,1),strong_ext10(2:end,2),'LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,4)
hold all 
hold all 
grid on
title('Strong Run-times 10-11')
xlabel('Cores');
ylabel('T_p (seconds');
xlim([1,48]);
xticks(cores);


errorbar(cores,strong_int11(:,1),strong_int11(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,strong_ext11(2:end,1),strong_ext11(2:end,2),'LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');





figure(3)

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
subplot(2,2,1)
hold all 
grid on
title('Strong Scaling 10-08')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');

errorbar(cores,int_scaling_strong08(:,1),int_scaling_strong08(:,2),'LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,ext_scaling_strong08(:,1),ext_scaling_strong08(:,2),'*','LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');


subplot(2,2,2)
hold all 
grid on
title('Strong Scaling 10-09')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');

errorbar(cores,int_scaling_strong09(:,1),int_scaling_strong09(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,ext_scaling_strong09(:,1),ext_scaling_strong09(:,2),'*','LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,3)
hold all 
grid on
title('Strong Scaling 10-10')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');

errorbar(cores,int_scaling_strong10(:,1),int_scaling_strong10(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,ext_scaling_strong10(:,1),ext_scaling_strong10(:,2),'*','LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');


subplot(2,2,4)
hold all 
grid on
title('Strong Scaling 10-11')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');

errorbar(cores,int_scaling_strong11(:,1),int_scaling_strong11(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,ext_scaling_strong11(:,1),ext_scaling_strong11(:,2),'*','LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');



figure(4)
set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
subplot(1,2,1)
hold all 
grid on
title('Strong Scaling, Internal Time (walltime)')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');

plot(cores,int_scaling_strong08(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,int_scaling_strong09(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,int_scaling_strong10(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,int_scaling_strong11(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
legend('N=10-08','N=10-09','N=10-10','N=10-11','Location','northwest','FontSize',20);

subplot(1,2,2)
hold all 
grid on
title('Strong Scaling, Ext Time (/usr/bin/time elapsed)')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');

plot(cores,ext_scaling_strong08(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,ext_scaling_strong09(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,ext_scaling_strong10(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,ext_scaling_strong11(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);

legend('N=10-08','N=10-09','N=10-10','N=10-11','Location','northwest','FontSize',20);


figure(5) 

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
hold all 
grid on
title('Strong Scaling, Internal Time (walltime)')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');

plot(cores,int_scaling_strong08(:,1),'k','HandleVisibility','off');
plot(cores,int_scaling_strong09(:,1),'k','HandleVisibility','off');
plot(cores,int_scaling_strong10(:,1),'k','HandleVisibility','off');
plot(cores,int_scaling_strong11(:,1),'k','HandleVisibility','off');










params = zeros(2,2,4);

[fit_int_8,fit_ext_8,params(:,:,1)] = fit_scaling(cores, int_scaling_strong08(:,1),ext_scaling_strong08(:,1) , serial_coef(2),  1e08);
[fit_int_9,fit_ext_9,params(:,:,2)] = fit_scaling(cores, int_scaling_strong09(:,1),ext_scaling_strong09(:,1) , serial_coef(2),1e09);
[fit_int_10,fit_ext_10,params(:,:,3)] = fit_scaling(cores, int_scaling_strong10(:,1),ext_scaling_strong10(:,1) , serial_coef(2),1e10);
[fit_int_11,fit_ext_11,params(:,:,4)] = fit_scaling(cores, int_scaling_strong11(:,1),ext_scaling_strong11(:,1) , serial_coef(2),1e11);

params




figure(6)

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
subplot(2,2,1)
hold all 
grid on
title('Strong Scaling 10-08')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');
errorbar(cores,int_scaling_strong08(:,1),int_scaling_strong08(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,ext_scaling_strong08(:,1),ext_scaling_strong08(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

plot(cores,fit_int_8,'LineWidth',1.8,'Color','b');
plot(cores,fit_ext_8,'LineWidth',1.8,'Color','r');
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,2)
hold all 
grid on
title('Strong Scaling 10-09')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');
errorbar(cores,int_scaling_strong09(:,1),int_scaling_strong09(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,ext_scaling_strong09(:,1),ext_scaling_strong08(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

plot(cores,fit_int_9,'LineWidth',1.8,'Color','b');
plot(cores,fit_ext_9,'LineWidth',1.8,'Color','r');
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,3)
hold all 
grid on
title('Strong Scaling 10-10')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');
errorbar(cores,int_scaling_strong10(:,1),int_scaling_strong10(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,ext_scaling_strong10(:,1),ext_scaling_strong10(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

plot(cores,fit_int_10,'LineWidth',1.8,'Color','b');
plot(cores,fit_ext_10,'LineWidth',1.8,'Color','r');
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,4)
hold all 
grid on
title('Strong Scaling 10-11')
xlabel('Cores');
ylabel('Speedup S');
xlim([1,48]);
xticks(cores);
plot(cores,cores,'k','HandleVisibility','off');
errorbar(cores,int_scaling_strong11(:,1),int_scaling_strong11(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,ext_scaling_strong11(:,1),ext_scaling_strong11(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

plot(cores,fit_int_11,'LineWidth',1.8,'Color','b');
plot(cores,fit_ext_11,'LineWidth',1.8,'Color','r');
legend('Walltime (internal)','Elapsed (external)','Location','northwest');





%parallel overhead


T_over_8 = parallel_overhead(cores,strong_ext08,strong_int08);
T_over_9 = parallel_overhead(cores,strong_ext09,strong_int09);
T_over_10 = parallel_overhead(cores,strong_ext10,strong_int10);
T_over_11 = parallel_overhead(cores,strong_ext11,strong_int11);


fit_over_8 = parallel_fit(cores,[serial_coef(2);params(:,1,1)], 1e08  );
fit_over_9 = parallel_fit(cores,[serial_coef(2);params(:,1,2)], 1e09  );
fit_over_10 = parallel_fit(cores,[serial_coef(2);params(:,1,3)], 1e10  );
fit_over_11 = parallel_fit(cores,[serial_coef(2);params(:,1,4)], 1e11  );


figure(7)

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
hold all 
grid on
title('Parallel Overhead fraction')
xlabel('Cores');
ylabel('Overhead');
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');
set(gca,'ColorOrderIndex',1)
errorbar(cores,T_over_8(:,1),T_over_8(:,2),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,T_over_9(:,1),T_over_9(:,2),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,T_over_10(:,1),T_over_10(:,2),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,T_over_11(:,1),T_over_11(:,2),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);

legend('N=10-08','N=10-09','N=10-10','N=10-11','Location','northwest','FontSize',20);


figure(8)

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
hold all 
grid on
title('Parallel Overhead fraction')
xlabel('Cores');
ylabel('Overhead');
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');
errorbar(cores,T_over_8(:,1),T_over_8(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,T_over_9(:,1),T_over_9(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,T_over_10(:,1),T_over_10(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,T_over_11(:,1),T_over_11(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

set(gca,'ColorOrderIndex',1)

plot(cores,fit_over_8,'LineWidth',2);
plot(cores,fit_over_9,'LineWidth',2);
plot(cores,fit_over_10,'LineWidth',2);
plot(cores,fit_over_11,'LineWidth',2);
legend('N=10-08','N=10-09','N=10-10','N=10-11','Location','northwest','FontSize',20);



figure(12)
hold all
set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);

title('Overhead Fitting coefficient k_1 vs Problem Size')
xlabel('Problem Size (N)');
ylabel('K_1 coefficient');


k1_f = @(c,x) transpose(c(1) + c(2).*x);
k1_coef0 = [0.001,1e-08];

k1_coef = lsqcurvefit(k1_f, k1_coef0, runs,...
    [params(1,1,1),params(1,1,2),params(1,1,3),params(1,1,4)],[0,0],[inf,inf])

k1_fit = k1_f(k1_coef,N_vec); 


errorbar(runs,...
    [params(1,1,1),params(1,1,2),params(1,1,3),params(1,1,4)],...
    [params(1,2,1),params(1,2,2),params(1,2,3),params(1,2,4)],...
    '*','LineStyle','none','Color','k','LineWidth',1.5,'MarkerSize',3);

%set(gca,'xscale','log');
%set(gca,'yscale','log');


plot(N_vec,k1_fit,'r','LineStyle',':','LineWidth',1.5);

axes('Position',[.15 .6 .3 .3])
box on
hold all
xlim([0.6e08,0.15e+10]);

errorbar(runs,...
    [params(1,1,1),params(1,1,2),params(1,1,3),params(1,1,4)],...
    [params(1,2,1),params(1,2,2),params(1,2,3),params(1,2,4)],...
    '*','LineStyle','none','Color','k','LineWidth',1.5,'MarkerSize',3);
plot(N_vec,k1_fit,'r','LineStyle',':','LineWidth',1.5);






ext_scaling_weak08 = scaling_err(serial_time_8,weak_ext08(:,:));
ext_scaling_weak09 = scaling_err(serial_time_9,weak_ext09(:,:));
ext_scaling_weak10 = scaling_err(serial_time_10,weak_ext10(:,:));
ext_scaling_weak11 = scaling_err(serial_time_11,weak_ext11(:,:));

int_scaling_weak08 = scaling_err(serial_time_8,weak_int08(:,:));
int_scaling_weak09 = scaling_err(serial_time_9,weak_int09(:,:,:));
int_scaling_weak10 = scaling_err(serial_time_10,weak_int10(:,:));
int_scaling_weak11 = scaling_err(serial_time_11,weak_int11(:,:));


% ext_scaling_weak08 = scaling_err(weak_ext08(1,:),weak_ext08(:,:));
% ext_scaling_weak09 = scaling_err(weak_ext09(1,:),weak_ext09(:,:));
% ext_scaling_weak10 = scaling_err(weak_ext10(1,:),weak_ext10(:,:));
% ext_scaling_weak11 = scaling_err(weak_ext11(1,:),weak_ext11(:,:));
% 
% int_scaling_weak08 = scaling_err(weak_int08(1,:),weak_int08(:,:));
% int_scaling_weak09 = scaling_err(weak_int09(1,:),weak_int09(:,:));
% int_scaling_weak10 = scaling_err(weak_int10(1,:),weak_int10(:,:));
% int_scaling_weak11 = scaling_err(weak_int11(1,:),weak_int11(:,:));





figure(33)

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
subplot(2,2,1)
hold all 
grid on
title('Weak Run-times 10-08')
xlabel('Cores');
ylabel('T_p (seconds');
xlim([1,48]);
xticks(cores);


errorbar(cores,weak_int08(:,1),weak_int08(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,weak_ext08(:,1),weak_ext08(:,2),'LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');


subplot(2,2,2)
hold all 
grid on
title('Weak Run-times 10-09')
xlabel('Cores');
ylabel('T_p (seconds');
xlim([1,48]);
xticks(cores);


errorbar(cores,weak_int09(:,1),weak_int09(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,weak_ext09(:,1),weak_ext09(:,2),'LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,3)
hold all 
hold all 
grid on
title('Weak Run-times 10-10')
xlabel('Cores');
ylabel('T_p (seconds');
xlim([1,48]);
xticks(cores);


errorbar(cores,weak_int10(:,1),weak_int10(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,weak_ext10(:,1),weak_ext10(:,2),'LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,4)
hold all 
hold all 
grid on
title('Weak Run-times 10-11')
xlabel('Cores');
ylabel('T_p (seconds');
xlim([1,48]);
xticks(cores);


errorbar(cores11,weak_int11(:,1),weak_int11(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores11,weak_ext11(:,1),weak_ext11(:,2),'LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');



figure(44)

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
subplot(2,2,1)
hold all 
grid on
title('Weak Scaling 10-08')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
yticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');

errorbar(cores,int_scaling_weak08(:,1),int_scaling_weak08(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,ext_scaling_weak08(:,1),ext_scaling_weak08(:,2),'*','LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');



subplot(2,2,2)
hold all 
grid on
title('Weak Scaling 10-09')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
yticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');


errorbar(cores,int_scaling_weak09(:,1),int_scaling_weak09(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,ext_scaling_weak09(:,1),ext_scaling_weak09(:,2),'*','LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,3)
hold all 
grid on
title('Weak Scaling 10-10')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
yticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');


errorbar(cores,int_scaling_weak10(:,1),int_scaling_weak10(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores,ext_scaling_weak10(:,1),ext_scaling_weak10(:,2),'*','LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');


subplot(2,2,4)
hold all 
grid on
title('Weak Scaling 10-11')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
yticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');


cores11 = [1,12,24,36,48];
errorbar(cores11,int_scaling_weak11(:,1),int_scaling_weak11(:,2),'*','LineStyle',':','Color','b','LineWidth',1.5,'MarkerSize',3);
errorbar(cores11,ext_scaling_weak11(:,1),ext_scaling_weak11(:,2),'*','LineStyle',':','Color','r','LineWidth',1.5,'MarkerSize',3);
legend('Walltime (internal)','Elapsed (external)','Location','northwest');




figure(55)
set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
subplot(1,2,1)
hold all 
grid on
title('Weak Scaling, Internal Time (walltime)')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
yticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');

plot(cores,int_scaling_weak08(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,int_scaling_weak08(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,int_scaling_weak10(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores11,int_scaling_weak11(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
legend('N=10-08','N=10-09','N=10-10','N=10-11','Location','northwest','FontSize',20);

subplot(1,2,2)
hold all 
grid on
title('Weak Scaling, Ext Time (/usr/bin/time elapsed)')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
yticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');

plot(cores,ext_scaling_weak08(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,ext_scaling_weak09(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores,ext_scaling_weak10(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);
plot(cores11,ext_scaling_weak11(:,1),'*','LineStyle','-','LineWidth',1.5,'MarkerSize',3);

legend('N=10-08','N=10-09','N=10-10','N=10-11','Location','northwest','FontSize',20);






fit_weak_int_08 = weak_scaling1(cores,[serial_coef(2);params(:,1,1)], 1e08  );
fit_weak_int_09 = weak_scaling1(cores,[serial_coef(2);params(:,1,2)], 1e09  );
fit_weak_int_10 = weak_scaling1(cores,[serial_coef(2);params(:,1,3)], 1e10  );
fit_weak_int_11 = weak_scaling1(cores,[serial_coef(2);params(:,1,4)], 1e11  );

fit_weak_ext_08 = weak_scaling2(cores,[serial_coef(2);params(:,1,1)], 1e08  );
fit_weak_ext_09 = weak_scaling2(cores,[serial_coef(2);params(:,1,2)], 1e09  );
fit_weak_ext_10 = weak_scaling2(cores,[serial_coef(2);params(:,1,3)], 1e10  );
fit_weak_ext_11 = weak_scaling2(cores,[serial_coef(2);params(:,1,4)], 1e11  );





figure(9)

set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
subplot(2,2,1)
hold all 
grid on
title('Strong Scaling 10-08')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');
errorbar(cores,int_scaling_weak08(:,1),int_scaling_weak08(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,ext_scaling_weak08(:,1),ext_scaling_weak08(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

plot(cores,fit_weak_int_08,'LineWidth',1.8,'Color','b');
plot(cores,fit_weak_ext_08,'LineWidth',1.8,'Color','r');
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,2)
hold all 
grid on
title('Strong Scaling 10-09')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');
errorbar(cores,int_scaling_weak09(:,1),int_scaling_weak09(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,ext_scaling_weak09(:,1),ext_scaling_weak09(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

plot(cores,fit_weak_int_09,'LineWidth',1.8,'Color','b');
plot(cores,fit_weak_ext_09,'LineWidth',1.8,'Color','r');
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,3)
hold all 
grid on
title('Strong Scaling 10-10')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');
errorbar(cores,int_scaling_weak10(:,1),int_scaling_weak10(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores,ext_scaling_weak10(:,1),ext_scaling_weak10(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

plot(cores,fit_weak_int_10,'LineWidth',1.8,'Color','b');
plot(cores,fit_weak_ext_10,'LineWidth',1.8,'Color','r');
legend('Walltime (internal)','Elapsed (external)','Location','northwest');

subplot(2,2,4)
hold all 
grid on
title('Strong Scaling 10-11')
xlabel('Cores');
ylabel('Efficiency E');
xlim([1,48]);
xticks(cores);
ylim([0.5,1.5]);
plot(cores,ones(length(cores),1),'k','HandleVisibility','off');
errorbar(cores11,int_scaling_weak11(:,1),int_scaling_weak11(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');
errorbar(cores11,ext_scaling_weak11(:,1),ext_scaling_weak11(:,2),'LineStyle','-','Color','k','LineWidth',1,'HandleVisibility','off');

plot(cores,fit_weak_int_11,'LineWidth',1.8,'Color','b');
plot(cores,fit_weak_ext_11,'LineWidth',1.8,'Color','r');
legend('Walltime (internal)','Elapsed (external)','Location','northwest');






function [scaling_out] = scaling_err(t_1,t_p)
    scaling_out = zeros(length(t_p(:,1)),2);
    
    scaling_out(:,1) = t_1(1)./t_p(:,1);
    err_1 = (t_1(2)/t_1(1))^2;
    for i=1:length(scaling_out(:,2))
        scaling_out(i,2) = sqrt( err_1 + (t_p(i,2)/t_p(i,1)).^2 )*abs(scaling_out(i,1));
    end
end


function [scaling_out] = th_scaling_noov(cores,params,N)

    s=params(1);
    p=params(2);
    
    cores_inv = 1./cores;
    
    
    scaling_out = (s + p*N)./( s + p*N.*cores_inv );
    %scaling_out=transpose(scaling_out);
    

end


function [scaling_out] = th_scaling_nocomm(cores,params,p,N)

    k1=params(1);
    
    cores_inv = 1./cores;
    
    %serial_t = s + p ;
    
    
   scaling_out = (N*p)./( p*N.*cores_inv + k1.*(cores - 1) );
   %scaling_out = (p)./( (p).*cores_inv + k1.*(1 -cores_inv));
    %scaling_out=transpose(scaling_out);
    

end

function [scaling_out] = th_scaling_comm(cores,params,p,N)

    k1=params(1);
    k2=params(2);
    
    cores_inv = 1./cores;
    
    %serial_t = s + p ;
    
    scaling_out = (p*N)./( N.*p.*cores_inv + k1.*(cores - 1)  + k2.*cores);
    %scaling_out = (p)./( (p).*cores_inv + k1.*(1 -cores_inv) + k2.*cores./N );
    %scaling_out=transpose(scaling_out);
    

end


function [out] = th_scaling_combined(cores,params,p,N)

    scaling_int = th_scaling_nocomm(cores,params,p,N);
    scaling_ext = th_scaling_comm(cores,params,p,N);
    
    out =  [scaling_int;scaling_ext] ;
    
end


function [scaling_int,scaling_ext,params] = fit_scaling(cores, t_internal,t_external ,p,N)

    x_data = transpose(cores);

    y_data =  [ t_internal ; t_external ];
    

    fun = @(params,x) th_scaling_combined(x,params,p,N);

    beta0 = [0.001,0.001];
    
    lb = [0,0];
    ub = [inf,inf];

    %[params,resnorm,residual,exitflag,output] = lsqcurvefit(fun, params0, x_data,y_data, lb, ub );
    
    [beta,R,J,sigma]  = nlinfit(x_data, y_data, fun, beta0);
    
    params=zeros(2:2);
    params(:,1) = beta;
    
    
    for i=1:2
      params(i,2) = sqrt(sigma(i,i));
    end
    
    params;
    
       

    scaling_int = th_scaling_nocomm(cores,params,p,N);
    scaling_ext = th_scaling_comm(cores,params,p,N);
    

    
end

function [T_over] = parallel_overhead(cores,external_t,internal_t)
    T_over = zeros(length(cores),2);
    t_serial =  external_t(1,1);
    t_serial_error =  external_t(1,2);
    
    
    external_t_error = external_t(2:end,2);
    external_t = external_t(2:end,1);
    
    internal_t_error = internal_t(2:end,2);
    internal_t = internal_t(2:end,1);
    
    for i = 1:length(cores)
        total_parallel_i =  external_t(i)*cores(i) ;
        total_parallel_i_err = external_t_error(i)*cores(i) ;
        
        T_over(i,1) = 1 - t_serial./total_parallel_i;
       
        T_over(i,2) = sqrt( (total_parallel_i_err/total_parallel_i)^2 + (t_serial_error/t_serial)^2 )*abs(T_over(i,1));
    end
   
    
%     for i = 1:length(cores)
%         total_parallel_i =  external_t(i)*cores(i) ;
%         total_parallel_i_err = external_t_error(i)*cores(i) ;
%         
%         T_over(i,1) = total_parallel_i - t_serial;
%        
%         T_over(i,2) = (total_parallel_i_err + t_serial_error)/T_over(i,1);
%         
%         T_over(i,1) = T_over(i,1) /total_parallel_i;
%         T_over(i,2) = sqrt( (total_parallel_i_err/total_parallel_i)^2 + (T_over(i,2))^2 )*abs(T_over(i,1));
%     end

end

function [over] = parallel_fit(cores,params,N)

    p=params(1);
    k1=params(2);
    k2=params(3);
    
   over =  1 - 1./( cores.*(cores.*(k1 + k2) - k1)/(p*N) + 1  );

end


function [scaling_out] = weak_scaling1(cores,params,N)

     p=params(1);
    k1=params(2);
    
    scaling_out = (p.*N )./( p.*N + k1.*(cores - 1).*cores );

end


function [scaling_out] = weak_scaling2(cores,params,N)

    params
     p=params(1);
    k1=params(2);
    k2=params(3);
    
    scaling_out = (p.*N )./( p.*N + k1.*(cores - 1).*cores + k2.*cores  );

end