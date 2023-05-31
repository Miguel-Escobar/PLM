using Plots, Graphs

function cargaInst(dirInst,instId)
   # Carga instancia de UFL dado el id de la instancia y
   # el directorio donde se encuentra.
   #

   # Carga matriz de coordenadas y costos de las localizaciones/facilities
   Fclt = readdlm(joinpath(dirInst,"fclt-" * instId * ".txt"), ',', Float64, '\n')
   # Carga matriz de coordenadas de clientes.
   Clnt = readdlm(joinpath(dirInst,"clnt-" * instId * ".txt"), ',', Float64, '\n')
   # Genera vector de costos y matriz de utilidad.
   n = size(Clnt,1)
   m = size(Fclt,1)
   costo = Fclt[:,3]

   util = zeros(n,m)
   for c in 1:n, l in 1:m
       util[c,l] = round.(10/((Fclt[l,1]-Clnt[c,1])^2+(Fclt[l,2]-Clnt[c,2])^2), digits=3)
   end

   return util, costo
end

function genArcos(x, util)
    # Entradas:
    #    x[]: vector indicatriz de facilities abiertas
    # Salida:
    #    arcos[]: vector de 2-tuplas [c,l] donde l es el indice
    #       de la localizacion a la que se conecta el cliente c.
    #
    arcos = []

    # Verificar si hay arcos y retornar si no.
    if sum(x) == 0
       return arcos
    end

    # Agregar arcos
    for c in 1:size(util,1)
    	append!( arcos, [[c, argmax( util[c,:] .* x)]] )
    end
    return arcos
end

function dibuja(Clnt, Fclt, arcos)
    # Entradas:
    #    Clnt: Matriz que describe clientes (coordenadas en el plano
    #	       cartesiano)
    #    Fclt: Matriz que describe clientes (coordenadas en el plano
    #	       cartesiano y costo de abrir localizaciones)
    #    arcos[]: vector de 2-tuplas [c,l] donde l es el indice
    #       de la localizacion a la que se conecta el cliente c.
    # Salida: Figura donde los clientes aparecen representados por círculos
    #         azules, las localizaciones como triángulos rojos, y los arcos
    #         corresponden a la asignación de clientes a localizaciones.
    #    
    n = size(Clnt,1)
    m = size(Fclt,1)
    scatter(cat(Fclt[:,1],Clnt[:,1],dims=1), cat(Fclt[:,2],Clnt[:,2],dims=1), markercolor = cat(fill(:red,m),fill(:blue,n),dims=1), markershape = cat(fill(:utriangle,m),fill(:circle,n),dims=1))
    for k in arcos
        plot!([Clnt[k[1],1], Fclt[k[2],1]],[Clnt[k[1],2],Fclt[k[2],2]], arrow = false, lc=:gray)
    end
    display(plot!(leg=false))
end
   