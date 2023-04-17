using JuMP, Gurobi

# Structure of a node in the tree (always copied). 
mutable struct BBNode
  # Tree structure. 
  parent::Union{Nothing, BBNode}
  children::Array{BBNode, 1}
  
  # Branching information. 
  variable::Int # Just an index of the variable for the branching.
                # Indexing starts at 1 (0: no branching; only for root node).
  variable_value::Int
  
  # Node's LP information.
  lp_value::Float64
  lp_feasibility::MOI.TerminationStatusCode   # :Optimal or something else. 
  solution::Array{Float64, 1}
  integer_feasible::Bool
end

mutable struct BBStatus
  upper_bound::Float64 # Computed from the LPs. 
  lower_bound::Float64 # Incumbent.
  incumbent::Array{Int, 1} # Rounded from floating-point solution. 
  
  root::BBNode
  # Optimisations: list of nodes to avoid going through all the tree,
  # sort nodes so the most interesting ones are picked up first, etc. 
end

