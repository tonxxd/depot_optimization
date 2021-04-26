
import numpy as np
import math
import copy
from mip import *

        

def parseTasks(TASKS, options):
    parsedTasks = copy.deepcopy(TASKS)
    slots = math.ceil(24*60/options["interval"])
    for task in parsedTasks:
        task["a"] =  math.ceil(task["a"] / options["interval"]) if math.ceil(task["a"] / options["interval"]) < slots else 0 # if depart after 23 index out of range
        task["d"] = math.floor(task["d"] / options["interval"])  
    return parsedTasks


def getPreconditioningAmount(t, TASKS, INTERVAL):
    amount = 0
    preconditioningSlots = math.ceil(60/INTERVAL) * 3
    for task in TASKS:
        if t >= task["d"]-preconditioningSlots and t < task["d"]-preconditioningSlots/3:
            amount = amount +.5
    return amount


def SchedulingModel(TASKS, upper_bound, userOptions = {}):
    

    options = userOptions
    
    
    # convert variables to correct interval
    TASKS = parseTasks(TASKS, options)
    
    m = Model(solver_name= GRB if options['solver'] == 'gurobi' else CBC)
    
    
    time = set(range(math.ceil(24*60/options["interval"])))
    tasks = set(range(len(TASKS)))

    
    PARAMS = {
        "TASKS": TASKS,
        "tasks": tasks,
        "upper_bound":upper_bound
    }
    
      ## PARAMS

    PARAMS["chargingRatePercent"] = [math.ceil(150*.9*.95*options["interval"]/60/TASKS[j]["batterySize"]*100) for j in tasks]
    PARAMS["last_slot"] = math.ceil(24*60/options["interval"]-1)
    
    PARAMS["preconditioning"] = [getPreconditioningAmount(t, TASKS, options["interval"]) for t in time]

    PARAMS["is_at_depot"] = [[True if ((t >= TASKS[j]['a'] or t < TASKS[j]['d']) and TASKS[j]['a']>TASKS[j]['d'] or ( t >= TASKS[j]['a'] and t < TASKS[j]['d'] and TASKS[j]['a']<TASKS[j]['d'])) else False for t in time] for j in tasks]

    PARAMS["cost"] = [1/(t+1) for t in time] if options["objective"] == 'cost' else [1 for t in time]
    
            
    ## VARS
    PARAMS["is_charging"] = [[m.add_var(var_type=BINARY, name=f'is_charging_{i}_{t}') for t in time] for i in tasks]
    PARAMS["Pc_grid"] = [m.add_var(var_type=CONTINUOUS,name=f'Pc_grid_{t}', lb=0, ub = len(TASKS)) for t in time]
    PARAMS["soc"] = [[m.add_var(var_type=CONTINUOUS,name=f'soc_{i}_{t}', lb=0, ub=100) for t in time] for i in tasks] 
    
    PARAMS["is_discharging"] = [[m.add_var(var_type=BINARY, name=f'is_discharging_{i}_{t}') for t in time] for i in tasks]
    
    #merge intervals together
    if not options["preemption"]:  
        PARAMS["is_charging_now_and_before"] = [[m.add_var(var_type=BINARY, name=f'is_charging_now_and_before_{i}_{t}') for t in time] for i in tasks]
    
    
    
    
    for t in time:
        
        m += PARAMS["Pc_grid"][t] == xsum(PARAMS["is_charging"][j][t] for j in tasks) - (xsum(PARAMS["is_discharging"][j][t] for j in tasks) if options["v2g"] else 0) + (PARAMS["preconditioning"][t] if options["preconditioning"] else 0), 'power balance'
        
        m+= PARAMS["Pc_grid"][t] <= upper_bound, 'max peak'
        
        for j in tasks:            

            if options["v2g"]:
                m += PARAMS["is_discharging"][j][t] + PARAMS["is_charging"][j][t] <= 1, 'discharge only if not charge'
                m += PARAMS["is_discharging"][j][t] <= PARAMS["is_at_depot"][j][t], 'discharge only if at depot'

            m += PARAMS["is_charging"][j][t] <= PARAMS["is_at_depot"][j][t], 'charge only if at depot'
            
            if not options["preemption"]:  
                m+=PARAMS["is_charging_now_and_before"][j][t] <= PARAMS["is_charging"][j][t], 'preemption constr 1'
                m+=PARAMS["is_charging_now_and_before"][j][t] <= PARAMS["is_charging"][j][PARAMS["last_slot"]] if t==0 else PARAMS["is_charging_now_and_before"][j][t]<=PARAMS["is_charging"][j][t-1]
                m+=PARAMS["is_charging_now_and_before"][j][t] >= PARAMS["is_charging"][j][PARAMS["last_slot"]] + PARAMS["is_charging"][j][t]-1 if t==0 else PARAMS["is_charging_now_and_before"][j][t] >= PARAMS["is_charging"][j][t-1] +PARAMS["is_charging"][j][t]-1
            
            if t == TASKS[j]['a']:
                m += PARAMS["soc"][j][t] == TASKS[j]['soc'], 'set soc at arrival'
              
            elif t == 0:
                if options["v2g"]:
                    m += PARAMS["soc"][j][t] == PARAMS["soc"][j][PARAMS["last_slot"]] + PARAMS["chargingRatePercent"][j] * (PARAMS["is_charging"][j][PARAMS["last_slot"]]-  PARAMS["is_discharging"][j][PARAMS["last_slot"]]), 'soc balance link'
                else:
                    m += PARAMS["soc"][j][t] == PARAMS["soc"][j][PARAMS["last_slot"]] + PARAMS["chargingRatePercent"][j] * (PARAMS["is_charging"][j][PARAMS["last_slot"]]), 'soc balance link'
            
            else:
                if options["v2g"]:
                    m += PARAMS["soc"][j][t] == PARAMS["soc"][j][t-1] + PARAMS["chargingRatePercent"][j] * (PARAMS["is_charging"][j][t-1] -  PARAMS["is_discharging"][j][t-1]),'soc balance'
                else:
                    m += PARAMS["soc"][j][t] == PARAMS["soc"][j][t-1] + PARAMS["chargingRatePercent"][j] * (PARAMS["is_charging"][j][t-1] ),'soc balance'

    for j in tasks:   
        m += PARAMS["soc"][j][int(TASKS[j]['d'])] >= 90, 'minimum charge at departure'
        
        if not options["preemption"]:  
            m += xsum(PARAMS["is_charging"][j][t] for t in time) - xsum(PARAMS["is_charging_now_and_before"][j][t] for t in time) <= options["limit_preemption"], 'limit_charging_intervals'



    #OBJECTIVE
    m.objective = xsum(PARAMS["cost"][t] * PARAMS["Pc_grid"][t] for t in time) 
    
    
    return m, PARAMS


def toOneMinuteInterval(arr, linear = False):
    finalArray = np.zeros(60*24)
    factor = int((60 * 24) / len(arr))
    for n in range(0, len(arr)):
        for i in range(n*factor, factor*(1+n)):
            base = i - n*factor
            finalArray[i] = arr[n] + base*(float(arr[n+1 if n+1 < len(arr) else 0] - arr[n])/float(factor)) if linear else arr[n]
    return list(finalArray);


def parse_results(modelVars):
    data = {
        "Pc_grid": list(map(lambda x: x.x,modelVars["Pc_grid"])), 
        "preconditioning":modelVars["preconditioning"], 
    }
    
    schedule = []
    for j in modelVars["tasks"]:        
        soc = list(map(lambda x: x.x, modelVars["soc"][j]))
        schedule.append({
            'bus':  modelVars["TASKS"][j]["name"],
            'arrival': modelVars["TASKS"][j]["a"],
            'departure': modelVars["TASKS"][j]["d"],
            "is_at_depot": list(map(lambda x: x, modelVars["is_at_depot"][j])),
            'is_charging': list(map(lambda x: x.x, modelVars["is_charging"][j])),
            'is_discharging': list(map(lambda x: x.x, modelVars["is_discharging"][j])),
            "soc":  list(map(lambda x: x, toOneMinuteInterval(soc, True))),
            "real_soc": soc
        })
    data["schedule"] = schedule
    return data

def get_min(TASKS,options={}, solver='cbc'):
    
    heuristic_starting_percent = 1
    current_peak = max(math.ceil(len(TASKS)*heuristic_starting_percent),1) 
    print(current_peak)
    if options["preconditioning"]:
        current_peak+=math.ceil(len(TASKS)*.5)
    print(current_peak)
    old_peak = current_peak
    
    lower_bound = 0
    iteration = 0
    optimalModel = False
    optimalPeak = 0
    optimalVars ={}
    optimalResults = False
    
    tried = []
    
    
    while(iteration == 0 or (current_peak not in tried and abs(old_peak-current_peak)>=1 and iteration < 10 and current_peak>=1)):
        
        print("TRYING",current_peak, solver)
        tried.append(current_peak)
    
        model, PARAMS = SchedulingModel(TASKS, upper_bound=current_peak, userOptions=options)
        model.emphasis=1
       
        results = model.optimize(max_nodes_same_incumbent = 1, max_solutions=1, max_seconds_same_incumbent=1)
        
        print("TRIED",current_peak, results)
    
        if(results == OptimizationStatus.OPTIMAL or results == OptimizationStatus.FEASIBLE):
            optimalModel = model
            optimalPeak = current_peak
            optimalVars = parse_results(PARAMS)
            optimalResults = results
            old_peak = current_peak
            current_peak = (lower_bound + (current_peak-lower_bound)/2) - (lower_bound + (current_peak-lower_bound)/2)%1
        else:
            lower_bound = current_peak
            current_peak = max((current_peak + (old_peak-current_peak)/2) - (current_peak + (old_peak-current_peak)/2)%1,1)
            if(lower_bound == current_peak):
                break

        iteration+=1
    
    
    
    if not optimalResults:
        raise Exception("Model is infeasible") 

    return optimalPeak, optimalModel, optimalResults, optimalVars


def get_power_profile(modelsVars, options):
    y = False
    for modelVars in modelsVars:
        y_model = modelVars["Pc_grid"]
        y_model = np.array(y_model)
        y = y_model if isinstance(y, bool) else np.sum([y, y_model], axis=0)
    return y




    
    
