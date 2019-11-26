function result = GAMultiObj(sysIn,is_pitch,is_first,first_td,settings,io)
%This function is defined to run the ga algorithm
sys = sysIn;

pop = 100;
if settings.Pop <= 2*pop
    pop = settings.Pop;
end

options = gaoptimset('Generations',settings.Gen,'PopulationSize',settings.Pop, 'Display', 'diagnose','UseParallel',true,'MigrationDirection','both','CrossoverFcn',@crossoverheuristic,'ParetoFraction',0.5,'PlotFcn',@gaplotpareto);

if is_pitch
    if settings.isDerEnable
        if is_first
            options = gaoptimset('Generations',settings.Gen,'PopulationSize',settings.Pop, 'Display', 'diagnose','UseParallel',true,'MigrationDirection','both','CrossoverFcn',@crossoverheuristic,'ParetoFraction',0.2,'PlotFcn',@gaplotpareto);
            lb = [settings.Kpmin;settings.Timin;settings.Kdmin;settings.Td];
            ub = [settings.Kpmax;settings.Timax;settings.Kdmax;2]; 
            [x,f,~] = gamultiobj(@(x)obj_p(x,sys,io),4,[],[],[],[],lb,ub,@(x)nonl_p(x,sys,io,settings),options);
        else
            options = gaoptimset('Generations',settings.Gen,'PopulationSize',pop, 'Display', 'diagnose','UseParallel',true,'MigrationDirection','both','CrossoverFcn',@crossoverheuristic,'ParetoFraction',0.5,'PlotFcn',@gaplotpareto);
            lb = [settings.Kpmin;settings.Timin;settings.Kdmin;first_td];
            ub = [settings.Kpmax;settings.Timax;settings.Kdmax;first_td]; 
            [x,f,~] = gamultiobj(@(x)obj_p(x,sys,io),4,[],[],[],[],lb,ub,@(x)nonl_p(x,sys,io,settings),options);
        end
    else
        lb = [settings.Kpmin;settings.Timin;settings.Kdmin;0];
        ub = [settings.Kpmax;settings.Timax;settings.Kdmax;0];
        
        if is_first
            options = gaoptimset('Generations',settings.Gen,'PopulationSize',settings.Pop, 'Display', 'diagnose','UseParallel',true,'MigrationDirection','both','CrossoverFcn',@crossoverheuristic,'ParetoFraction',0.2,'PlotFcn',@gaplotpareto);
        else
            options = gaoptimset('Generations',settings.Gen,'PopulationSize',pop, 'Display', 'diagnose','UseParallel',true,'MigrationDirection','both','CrossoverFcn',@crossoverheuristic,'ParetoFraction',0.5,'PlotFcn',@gaplotpareto);
        end
        
        [x,f,~] = gamultiobj(@(x)obj_p(x,sys,io),4,[],[],[],[],lb,ub,@(x)nonl_p(x,sys,io,settings),options);
    end
else
    options = gaoptimset('Generations',settings.Gen,'PopulationSize',pop, 'Display', 'diagnose','UseParallel',true,'MigrationDirection','both','CrossoverFcn',@crossoverheuristic,'ParetoFraction',0.5,'PlotFcn',@gaplotpareto);
    lb = [settings.Qkpmin;settings.Qtimin];
    ub = [settings.Qkpmax;settings.Qtimax];
    [x,f,~] = gamultiobj(@(x)obj_q(x,sys,io),2,[],[],[],[],lb,ub,@(x)nonl_q(x,sys,io,settings),options);
end

result = [x,f];

end
    




