clc;close all;

cores=[1,4,8,12,16,20,24,28,32,36,40,44,48];



%parallel overhead


T_over_8 = parallel_overhead2(cores,elapsed_08(1:end,1),internal_08(1:end,1));
T_over_9 = parallel_overhead2(cores,elapsed_09(1:end,1),internal_09(1:end,1));
T_over_10 = parallel_overhead2(cores,elapsed_10(1:end,1),internal_10(1:end,1));
T_over_11 = parallel_overhead2(cores,elapsed_11(1:end,1),internal_11(1:end,1));


figure(7)
hold all 

plot(cores,T_over_8);
plot(cores,T_over_9);
plot(cores,T_over_10);
plot(cores,T_over_11);



function [T_over] = parallel_overhead(cores,external_t,internal_t)
    T_over = zeros(length(cores),1);
    
    T_mpi = external_t(2) - external_t(1);
    
    external_t = external_t(2:end);
    
    for i = 1:length(cores)
        T_over(i) = T_mpi + cores(i)*internal_t(i) - internal_t(1);
        T_over(i) = T_over(i)./external_t(i);
    end
    
    internal_t
    T_over

end

function [T_over] = parallel_overhead2(cores,external_t,internal_t)
    T_over = zeros(length(cores),1);
    t_serial =  external_t(1);
    
     external_t = external_t(2:end);
    
    T_mpi = external_t(1) - t_serial;
    
   
    
    for i = 1:length(cores)
        total_parallel_i =  external_t(i)*cores(i) ;
        
        T_over(i) = total_parallel_i - t_serial;
        T_over(i) = T_over(i) /total_parallel_i;
    end
    
 

end