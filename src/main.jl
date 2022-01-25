
action = Args.get(:action)

action == :dual           && @ms include("dual.jl")
action == :plans_coupants && @ms include("plans_coupants.jl")
action == :branch_cut     && @ms include("branch_cut.jl")
action == :heuristic      && @ms include("heuristic.jl")


function main()

    time1 = time()

    action = Args.get(:action)
    println("action = $action")

    if action == :dual      

    
    elseif action == :plans_coupants 
    
    
    elseif action == :branch_cut     
    
    
    else if action == :heuristic 

    elseif action == :none
        println("Aucune action indiquée")
        println(Args.get_syntaxe())
    else
        println("Erreur : action $(action) non implémentée (dans main.jl)")
        println(Args.get_syntaxe())
        exit(1)
    end

    time2 = time()
    sec = round((time() - time1), digits = 3) # on veut limiter la précision à la ms
    ln1("Durée totale du main 1000*(time2-time1) : $(sec)s")
    ln1("main() END")
end
        