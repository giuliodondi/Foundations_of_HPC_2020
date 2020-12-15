clear;clc;close all;

format long
figure_size=1024;

t_c = 2*10^(-9);
t_r = 1*10^(-4);
t_comm = 1*10^(-6);

cores = 1:1:100;





%N=[100000,500000,1000000,5000000,10000000,50000000,100000000];
N = [20000,100000,200000,1000000,20000000];

figure(1)
hold all
grid on
title('Speedup without comm.overhead')
set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
set(gca, 'OuterPosition', [0,0,1,1])
xlabel('Cores');
ylabel('Speedup S');
set(gca,'yscale','log')
plot(cores,cores,'k','HandleVisibility','off');
speedup = zeros(length(cores),1);

for i=1:length(N)

   s_t = serial_time(N(i),t_c,t_r);
   for j=1:length(cores)
       speedup(j) = s_t/parallel_time(N(i),cores(j),t_c,t_r,0);
   end
   
   plot(cores,speedup,'LineWidth',1.5);
    
end

legend('N=20000','N=100000','N=200000','N=1000000','N=20000000','Location','northwest');

figure(2)
hold all
grid on
title('Speedup with overhead')
set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
set(gca, 'OuterPosition', [0,0,1,1])
xlabel('Cores');
ylabel('Speedup S');
set(gca,'yscale','log')
plot(cores,cores,'k','HandleVisibility','off');
speedup = zeros(length(cores),1);



for i=1:length(N)

   s_t = serial_time(N(i),t_c,t_r);
   for j=1:length(cores)
       speedup(j) = s_t/parallel_time(N(i),cores(j),t_c,t_r,t_comm);
   end
   
   plot(cores,speedup,'LineWidth',1.5);
    
end

legend('N=20000','N=100000','N=200000','N=1000000','N=20000000','Location','northwest');


figure(3)
hold all
grid on
title('Speedup, alternative algorithm')
set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
xlabel('Cores');
ylabel('Speedup S');
set(gca,'yscale','log')
plot(cores,cores,'k','HandleVisibility','off');
speedup = zeros(length(cores),1);


for i=1:length(N)

   s_t = serial_time(N(i),t_c,t_r);
   for j=1:length(cores)
       speedup(j) = s_t/parallel_time(N(i),cores(j),t_c,t_r,t_comm);
   end
   
   [maxx, idxx] = max(speedup);
  
   
   disp(['for N= ', num2str(N(i)), ' max is at ', num2str(cores(idxx)) ,' cores.']);
   
  plot(cores,speedup,'--','Color','k','HandleVisibility','off');
    
end
set(gca,'ColorOrderIndex',1)

for i=1:length(N)

   s_t = serial_time(N(i),t_c,t_r);

   for j=1:length(cores)
       speedup(j) = s_t/parallel_time_improved(N(i),cores(j),t_c,t_r,t_comm);
   end
   
     [maxx, idxx] = max(speedup);
   disp(['for N= ', num2str(N(i)), ' max is at ', num2str(cores(idxx)) ,' cores.']);
   
   plot(cores,speedup,'LineWidth',1.5);
end
legend('N=20000','N=100000','N=200000','N=1000000','N=20000000','Location','northwest');


figure(30)
hold all
grid on
title('Speedup, alternative algorithm')
set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
xlabel('Cores');
ylabel('Speedup S');
set(gca,'yscale','log')
plot(cores,cores,'k','HandleVisibility','off');
speedup = zeros(length(cores),1);


for i=1:length(N)

   s_t = serial_time(N(i),t_c,t_r);
   for j=1:length(cores)
       speedup(j) = s_t/parallel_time(N(i),cores(j),t_c,t_r,t_comm);
   end
   
   [maxx, idxx] = max(speedup);
  
   
   disp(['for N= ', num2str(N(i)), ' max is at ', num2str(cores(idxx)) ,' cores.']);
   
  plot(cores,speedup,'--','Color','k','HandleVisibility','off');
    
end

for i=1:length(N)

   s_t = serial_time(N(i),t_c,t_r);

   for j=1:length(cores)
       speedup(j) = s_t/parallel_time_improved2(N(i),cores(j),t_c,t_r,t_comm);
   end
   
     [maxx, idxx] = max(speedup);
   disp(['for N= ', num2str(N(i)), ' max is at ', num2str(cores(idxx)) ,' cores.']);
   
   plot(cores,speedup,'LineWidth',1.5);
end
legend('N=20000','N=100000','N=200000','N=1000000','N=20000000','Location','northwest');




figure(4)
hold all
grid on
title('Speedup differences')
set(gcf,'units','pixel','position',[0,0,figure_size,figure_size]);
xlabel('Cores');
ylabel('Speedup difference');
speedup = zeros(length(cores),1);

for i=1:length(N)

   s_t = serial_time(N(i),t_c,t_r);

   for j=1:length(cores)
       speedup(j) = s_t*(1./parallel_time_improved(N(i),cores(j),t_c,t_r,t_comm) - 1./parallel_time(N(i),cores(j),t_c,t_r,t_comm) );
   end
   
   
   plot(cores,speedup,'LineWidth',1.5);
end
legend('N=20000','N=100000','N=200000','N=1000000','N=20000000','Location','northwest');


function [time] = parallel_time(N,p, t_c,t_r,t_comm)

task_frac = ceil(N/p-1);

time = t_r + t_comm*2*(p-1) + t_c*(p - 2 + task_frac);

end


function [time] = serial_time(N, t_c,t_r)

time = t_r + t_c*(N - 1);

end

function [time] = parallel_time_improved(N,p, t_c,t_r,t_comm)

task_frac = ceil(N/p-1);

time = t_r + t_comm*(p-1) + t_c*(task_frac - 1) + (t_c + t_comm)*(log(p)./log(2));

end


function [time] = parallel_time_improved2(N,p, t_c,t_r,t_comm)

task_frac = ceil(N/p-1);

time = t_r + t_comm*(p-1) + t_c*(task_frac - 1) + (t_c + t_comm)*(log(p)./log(2) + mod(str2double(dec2base(p,2)),9) - 1);

end
