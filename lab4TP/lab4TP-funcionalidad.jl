using Plots, Graphs

function cargaInst(dirInst,instId)
   # Carga instancia de UFL dado el id "instId" de la instancia y
   # el directorio "dirInst" donde se encuentra.
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

   return util, costo,  Fclt, Clnt
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
    	maxUtil = -Inf
	maxL = 0
    	for l in 1:size(util,2)
	    if x[l] == 1 && util[c,l] > maxUtil
	       maxUtil = util[c,l]
	       maxL = l
	    end
	end
	append!( arcos, [[c, maxL]] )
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
    scatter!(cat(Fclt[:,1],Clnt[:,1],dims=1), cat(Fclt[:,2],Clnt[:,2],dims=1), markercolor = cat(fill(:red,m),fill(:blue,n),dims=1), markershape = cat(fill(:utriangle,m),fill(:circle,n),dims=1))
    display(plot!(leg=false,size=(900,600),thickness_scaling=0.5))
end

function initMu(util)
    # Entrada:
    #     util: Matriz de utilidad de UFL,
    # Salida: Vector de valores iniciales para la relajación lagrangiana LR(mu)
    #
    cardC = size(util,1);

    mu = zeros(cardC);
    for c in 1:cardC
        mu[c] = (reverse(sort(util[c,:])))[2];
    end

    return mu
end

function cotaInf(mu, util, costo)
    # 
    s = sum( max.(zeros(size(util)), util - mu*ones(size(util,2))'), dims=1)
    x = zeros(size(util,2))
        
    for l in 1:size(util,2)
        if s[l] > costo[l] - 0.1
            x[l] = 1
        end
    end

    return x
end

function zUFL(x, util, costo)
    # Entradas:
    #   x: Vector indicatriz de localizaciones abiertas
    #   util: Matriz de utilidad de UFL,
    #   costo: Vector de costo de habilitar localizaciones
    # Salida: Vector de valores iniciales para la relajación lagrangiana LR(mu)
    #
    cardC = size(util,1)
    cardL = size(util,2)

    maxU = -Inf*ones(cardC);
    for c in 1:cardC
        for l in 1:cardL
            if ( x[l] > 0.99 ) && ( util[c,l] > maxU[c] ) 
                maxU[c] = util[c,l]
            end
         end
    end

    val = sum( maxU )
    for l in 1:cardL
    	val = val - x[l] * costo[l]
    end

    return val
end

function zLR(mu, util, costo)
    #
    #
    #
    cardC = size(util,1)
    cardL = size(util,2)

    foo = zeros(cardL)
    for l in 1:cardL, c in 1:cardC
        foo[l] = foo[l] + max(0.0, util[c,l] - mu[c])
    end 
    
    obj_zLR = sum( mu )
    for l in 1:cardL
        obj_zLR = obj_zLR + max(0.0, foo[l] - costo[l])
    end
    
    return obj_zLR
end