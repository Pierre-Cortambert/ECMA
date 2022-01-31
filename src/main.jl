using JuMP
using CPLEX
using ArgParse

include("read_write.jl")
include("static.jl")
include("dual.jl")
include("plans_coupants.jl")
include("plan_coupant_knackpack.jl")
include("branch_cut.jl")
include("heuristic.jl")


function parse_commandline()
    settings = ArgParseSettings()
    @add_arg_table settings begin
        "--algo"
            help = "Algorithme à choisir (parmi : s, bc, d, h, pc, pck)"
            arg_type = String
            default = "none"
        "--file"
            help = "Nom da l'instance à considérer"
            arg_type = String
            default = "20_USA-road-d.BAY.gr"
        "--time"
            help = "Temps limite de recherche de solution optimale"
            arg_type = Float64
            default = Inf # pas de limitation
        "--pr" #critère d'arret correspondant au pourcentage de la solution robuste trouvé par rapport à la solution statique  : non implémentée encore
            default = 0
    end
    return parse_args(ARGS, settings)
end


function main()
    parsed_args = parse_commandline()
    println("../instances/"*parsed_args["file"])
    (n,s,t,S,d1,d2,p,ph,A,d,D)= get_data("../instances/"*parsed_args["file"])
    time1 = time()

    algo = parsed_args["algo"]
    if algo == "s"
        isOptimal, traj, sol = Static(n,s,t,S,d1,d2,p,ph,A,d,D)
        cpt = 0

    elseif algo == "d"      
        isOptimal, traj, sol = Dual(n,s,t,S,d1,d2,p,ph,A,d,D) 
        cpt = 0
    
    elseif algo == "pc"
        isOptimal, traj, sol, cpt = plan_coupant(n,s,t,S,d1,d2,p,ph,A,d,D,parsed_args["time"])   
    
    elseif algo == "pck"
        isOptimal, traj, sol, cpt = plan_coupant_kp(n,s,t,S,d1,d2,p,ph,A,d,D,parsed_args["time"])  

    elseif algo == "bc"    
        
    elseif algo == :"h"

    elseif action == "none"
        println("Aucune action indiquée")
        println(Args.get_syntaxe())
    else
        println("Erreur : action $(action) non implémentée (dans main.jl)")
        return
    end
    time2 = time()
    sec = round((time() - time1), digits = 3) # on veut limiter la précision à la ms
    println("Durée totale du main 1000*(time2-time1) : $(sec)s")
    write_sol(algo,parsed_args["file"],isOptimal, traj, sol, cpt , sec)
    
end
        
main()