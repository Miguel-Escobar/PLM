n = 8

fil = [];
col = [];
cap = [];
for k in 1:4
  d = (k-1)*n
  for i in 1:n
    for j in i+1:n
      if rand() < 2/3
        push!( fil, i+d )
	push!( col, j+d )
        push!( cap, round(rand(Float64),digits=2) )
      end
    end
  end
end

for k in 1:4
  for m in 1:4
    for t in 1:2
      i = rand(1:n)
      j = rand(1:n)
      if ( i+(k-1)*n < j+(m-1)*n )
        r = round(rand(Float64),digits=2)
	push!( fil, i+(k-1)*n)
	push!( col, j+(m-1)*n)
	push!( cap, r )
      end
      if ( i+(k-1)*n > j+(m-1)*n )
        r = round(rand(Float64),digits=2)
	push!( fil, j+(m-1)*n)
	push!( col, i+(k-1)*n)
	push!( cap, r )
      end
    end
  end
end

using Printf;
for i in 1:length(fil)
  @printf( "%d ", fil[i] )
end
@printf( "\n" )
for i in 1:length(fil)
  @printf( "%d ", col[i] )
end
@printf( "\n" )
for i in 1:length(fil)
  @printf( "%.2f ", cap[i] )
end
@printf( "\n" )



###############################
using Graphs, GraphsFlows
flow_graph = Graphs.DiGraph(4*n) # Create a flow graph
for i in 1:4*n
  for j in i+1:4*n
    if U[i,j] != 0
        Graphs.add_edge!(flow_graph, i, j)
        Graphs.add_edge!(flow_graph, j, i)
    end
  end
end

t = 3
S, noS, f = maximum_flow(flow_graph, s, t, U, algorithm=BoykovKolmogorovAlgorithm());