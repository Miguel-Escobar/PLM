# Lab 2, PLM 2023
#
using JuMP, Gurobi

###################################
# INCLUDE DESCRITION OF STRUCTURES
include("lab2TI-Structs.jl");

###################################
# GLOBAL VARIABLES
W = [];
const GUROBI_ENV = Gurobi.Env()      #Esta referencia  nos servirá
                                     #para usar solo un ambiente de Gurobi.

###################################
# FUNCTIONS 
function solver()
  mdl = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0, "Presolve" => 0));
  return(mdl)
end

function get_lp()
  mdl = solver();

  wght = W[:,1];
  valr = W[:,2];
  n = length(valr);
    
  @variable(mdl, 0 <= x[1:n] <= 1);
  @constraint(mdl, sum(wght .* x) <= n*5);
  @objective(mdl, Max, sum(valr .* x) );
  return mdl, all_variables(mdl)
end

integer_precision = 1e-4; # Integers are represented as floating-point
                          # values: if they are close enough to an integer,
			  # they are considered integer.
gap_precision = 1e-3;

function solution_is_boolean(sol::Array{Float64,1})
  for s in sol
    if s > integer_precision && s < 1 - integer_precision
      return false
    end
  end
  return true
end

function gap(status::BBStatus)
  # Will return NaN if no integer-feasible solution is found. 
  return (status.upper_bound - status.lower_bound) / status.upper_bound
end

function is_pruned(node::BBNode, status::BBStatus; quiet::Bool = true)
  if node.lp_feasibility != MOI.OPTIMAL
    if ! quiet println("out Node pruned by infeasibility.") end
    return true, :Feasibility
  end
  
  if node.integer_feasible
    if ! quiet println("out Node pruned by optimality.") end
    return true, :Optimality
  end

  if node.lp_value < status.lower_bound
    if ! quiet printn("out Node pruned by bounds.") end
    return true, :Bounds
  end
  
  return false, :None
end

function get_next_node(status)
  # Find a node which has no children and is not pruned. 
  # Pruning is then used here, together with tree-exploration strategy! 
  
  function explore_node(node::BBNode)
    pruned, why_pruned = is_pruned(node, status; quiet=true)
    if pruned
      return nothing
    end
  
    if length(node.children) == 0
      #println("out self no children: " * string(node))
      return node
    end
    
    # Binary variables: two children or zero, no other case. 
    node_left = node.children[1]
    if node_left != nothing
      possible_left = explore_node(node_left)
    else
      possible_left = nothing
    end
    
    node_right = node.children[2]
    if node_right != nothing
      possible_right = explore_node(node_right)
    else
      possible_right = nothing
    end
    
    if possible_left != nothing
      return possible_left
    elseif possible_right != nothing
      return possible_right
    end
  end
  
  return explore_node(status.root)
end

function is_variable_branched_on(node::BBNode, var::Int)
  # Depends on the path followed from root to now. 
  if node.variable == var
    return true
  end 
  
  if node.variable == 0 # Root node. 
    return false
  end
  
  return is_variable_branched_on(node.parent, var)
end

function find_branching_variable(sol::Array{Float64,1}, current_node::BBNode)
  # Only called on non-integer solutions. 
  n_vars = length(sol)
  for s in 1:n_vars
    v = sol[s]
    if v > integer_precision && v < 1 - integer_precision && ! is_variable_branched_on(current_node, s)
      return s
    end
  end
  println(" -> Where shall I branch?")
end

function branch_lp(current_node::BBNode, new_variable::Int, new_value::Int)
  lp, vars = get_lp()
  
  # Retrieve the previous nodes' branching decisions. 
  node = current_node
  while true
    if node.variable == 0
      break
    end
    
    @constraint(lp, vars[node.variable] == node.variable_value)
    node = node.parent
  end
  
  @constraint(lp, vars[new_variable] == new_value)
  return lp, vars
end

function create_node(parent::BBNode, new_variable::Int, new_value::Int)
  lp, lp_vars = branch_lp(parent, new_variable, new_value)
  optimize!(lp);
  lp_feasibility = termination_status(lp);
  lp_value = objective_value(lp);
  lp_solution = map(x->value(x), lp_vars);
  int_feas = solution_is_boolean(lp_solution);
  return BBNode(parent, [], new_variable, new_value, lp_value, lp_feasibility, lp_solution, int_feas);
end

function round_solution(sol::Array{Float64, 1})
  return [(s < .5) ? 0 : 1 for s in sol]
end

function get_upper_bound(status::BBStatus)
  # As the tree is developed such that both children are always created at the same time, 
  # retrieve the leaves' LP values, take the minimum of them. 
  function get_leaves(node::BBNode)
    if length(node.children) == 0
      pruned, why = is_pruned(node, status)
      if node.integer_feasible || (pruned && why != :Bounds)
        return -Inf
      else
        return node.lp_value
      end
    else
      left = get_leaves(node.children[1])
      right = get_leaves(node.children[2])
      return max(left, right)   
    end
  end
    
  return get_leaves(status.root)
end

import Base.print

function print(tree::BBNode, status::BBStatus)
  function space(level::Int)
    return "  " ^ level
  end

  function sub_print(node::BBNode, level::Int)
    print("out ")
    print(space(level))
    
    if node.variable != 0
      print("Variable " * string(node.variable) * "   :   " * string(node.variable_value) * "   ;   " * string(node.lp_value))
    else
      print("Root node:               " * string(node.lp_value))
    end
    
    if node.integer_feasible
      print("  (I)")
    else
      print("     ")
    end
    
    pruned, why_pruned = is_pruned(node, status; quiet=true)
    if pruned
      print("  (P: " * string(why_pruned) * ")")
    else
      print("     ")
    end
    
    print("\n")
    
    for child in node.children
      sub_print(child, level + 1)
    end
  end
  
  sub_print(tree, 2)
end

function print(status::BBStatus)
  println("out Complete branch-and-bound status, including tree.")
  println("out Upper bound: " * string(status.upper_bound))
  println("out Lower bound: " * string(status.lower_bound))
  println("out Incumbent: " * string(status.incumbent))
  println("out Tree: ")
  print(status.root, status)
end

function bnb()
  root_lp, root_vars = get_lp()
  optimize!(root_lp)
  root_feasibility = termination_status(root_lp);
    # Return values:
    #   http://www.juliaopt.org/JuMP.jl/v0.19/solutions/#JuMP.termination_status
  if  termination_status(root_lp) != MOI.OPTIMAL
    if termination_status(root_lp) == MOI.INFEASIBLE
      println("out Root node could not be solved to optimality: Infeasible" )
    elseif termination_status(root_lp) == MOI.DUAL_INFEASIBLE
      println("out Root node could not be solved to optimality: Dual_Infeasible" )
    else
      println("out Root node could not be solved to optimality: Undef" )
    end
    return
  end
  
  root_solution = map(x->value(x), root_vars);
  root_obj = objective_value(root_lp);
  if solution_is_boolean(map(x->value(x),root_vars))
    println("out Solved at root node.")
    println("out ", map(x->value(x),root_vars))
    println("out Objective: " * string(root_obj))
    root = BBNode(nothing, [], 0, 0, root_obj, root_feasibility, root_solution, true)
    status = BBStatus(root_obj, root_obj, map(x->value(x),root_vars), root)
    return status
  end

  root = BBNode(nothing, [], 0, 0, root_obj, root_feasibility, root_solution, false)
  status = BBStatus(root_obj, -Inf, [], root)
  
  println("out [0] Root node solved. LP value: " * string(root_obj) * ".")
  
  # Node exploring: do both children at once, while exploring current node.
  # Not what is usually done, but easier (no list of nodes to explore). 
  n_nodes_explored = 0

  println( gap(status) )
  while gap(status) > gap_precision
    n_nodes_explored += 1
    println("out [" * string(n_nodes_explored) * "] Exploring one more node.")
    
    # Find a node to explore. There may be none. 
    node_to_explore = get_next_node(status)
    if node_to_explore == nothing
      println("out No more branches to explore in the tree.")
      break
    end

    # Decide the branching variable. 
    var_index = find_branching_variable(node_to_explore.solution, node_to_explore)
    println("out [" * string(n_nodes_explored) * "] Branching on variable " * string(var_index))
    
    # Create new nodes with branching on this *binary* variable (with all constraints from root to this new branch).
    left = create_node(node_to_explore, var_index, 0)
    right = create_node(node_to_explore, var_index, 1)

    # Update the parent with the new children. 
    node_to_explore.children = [left, right]
    
    # Update the bounds and, if necessary, update the incumbent. 
    if left.integer_feasible && left.lp_value > status.lower_bound
      println("out [" * string(n_nodes_explored) * "] Found a new incumbent. Value: " * string(left.lp_value))
      status.lower_bound = left.lp_value
      status.incumbent = round_solution(left.solution)
    end
    if right.integer_feasible && right.lp_value > status.lower_bound
      println("out [" * string(n_nodes_explored) * "] Found a new incumbent. Value: " * string(right.lp_value))
      status.lower_bound = right.lp_value
      status.incumbent = round_solution(right.solution)
    end
    status.upper_bound = get_upper_bound(status)
    
    # End iteration. 
    println("out [" * string(n_nodes_explored) * "] Iteration done. Gap: " * string(gap(status)) * ". Bound: [" * string(status.lower_bound) * "; " * string(status.upper_bound) * "].")
  end

  return status
end
