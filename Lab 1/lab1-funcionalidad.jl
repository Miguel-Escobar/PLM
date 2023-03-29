### Archivo lab1-funcionalidad.jl
### Version 2023, curso MA4702.
### Universidad de Chile.


#Funciones auxiliares

function dibuja(coordx,coordy,arcos)
# Recibe dos arreglos de N valores donde (coordx[i],coordy[i])
# son las coordenadas del punto i.
# Recibe además una matriz arcos de N x N, donde arcos[i,j] es el
# peso del arco [i,j]
# Dibuja los N puntos en el plano y dibuja los arcos con ancho
# de linea proporcional al peso.
#
N=size(coordx,1)
    scatter(coordx,coordy,txt=text.(1:N,10,:bottom))
    for k in findall(!iszero, arcos)
        plot!([coordx[k[1]],coordx[k[2]]],[coordy[k[1]],coordy[k[2]]],arrow = false,lc=:blue, linewidth = 2*arcos[k[1],k[2]])
    end
    display(plot!(leg=false))
end

function(lee_archivo(nombre_archivo))
# Recibe un archivo con las coordenadas de N puntos.
# Devuelve:
#   - No. de puntos N,
#   - coordenadas x,
#   - coordenadas y,
#   - cota en el grado para el punto
#
    archivo = open(nombre_archivo)
    Lineas = readlines(archivo)
    n = size(Lineas,1)
        coordx = []
        coordy = []
	cotab = []
    for i in 1:n
        x,y,b=split(Lineas[i])
        push!(coordx,parse(Float64,x))
        push!(coordy,parse(Float64,y))
        push!(cotab,parse(Float64,b))
    end
    return n,coordx,coordy,cotab
end;

function componente_conexa(A)
# Recibe una matriz de adyacencia A de un digrafo G, ignora la orientación
# de los arcos y retorna un arreglo cuyos elementos son arreglos de enteros
# (tantos elementos como componentes conexas tenga el digrafo G).
# Cada uno de los arreglos de enteros corresponde a los índices de los
# vértices de una componente conexa de G.
#
  eps = 0.00001
  g = SimpleGraph(A .>= 1-eps);
  cmpnt = weakly_connected_components(g);
  return cmpnt
end