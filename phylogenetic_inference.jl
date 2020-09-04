module TreeDefinition
mutable struct PhyloTree
    name::AbstractString
    dist::Float64
    children::Vector{PhyloTree}
    mother::Union{PhyloTree, Missing}
end

function PhyloTree(name::String=nm, dist=1.0)
    PhyloTree(name, dist, PhyloTree[], missing)
end

export PhyloTree
end

using .TreeDefinition

function add_child!(a::PhyloTree, b::PhyloTree)
    push!(a.children, b)
    b.mother = a
end

# Prints a horizontal representation of the tree with centered branches.
# Example:
#       /-Spanish
#    /-|
#   |   \-Italian
#   |
# --|         /-Dutch
#   |      /-|
#   |   /-|   \-German
#   |  |  |
#    \-|   \-English
#      |
#       \-Gothic
function print(tree::PhyloTree)
    for line in tree_to_lines(tree, "--")[1]
        println(line)
    end
end


# Helper function for `print`.
# Arguments:
# - `current`: A PhyloTree.
# - `prefix`: The two-character prefix of the current leaf node / branch / tree.
# Returns:
# - A vector of lines to be printed
# - The index of the left-most branch connector (relevant for recursive calls)
function tree_to_lines(current::PhyloTree, prefix::String)
    if isempty(current.children)
        # Reached a leaf node.
        # The vertical centre of a leaf node is the leaf node itself (index 1).
        return ([prefix * current.name], 1)
    end

    lines = []
    # 'Pipe' refers to the vertical lines that connect sibling branches.
    pipe_start = -1
    pipe_end = -1
    uppermost_child = current.children[1]
    lowermost_child = current.children[length(current.children)]
    for child in current.children
        # Upper branch / lower branch / middle branch (if tree is not binary)
        if child == uppermost_child
            branch_pfx = "/-"
        elseif child == lowermost_child
            branch_pfx = "\\-"
        else
            branch_pfx = "--"
        end

        child_branch = tree_to_lines(child, branch_pfx)
        # Skip all previous lines and get the indices of the lines below(/above)
        # the uppermost(/lowermost) child to mark the beginning/end of the
        # vertical sibling connection line.
        if child == uppermost_child
            pipe_start = length(lines) + child_branch[2] + 1
        elseif child == lowermost_child
            pipe_end = length(lines) + child_branch[2] - 1
        end

        append!(lines, child_branch[1])

        # Add a line of vertical space between each pair of node names.
        if child != lowermost_child
            push!(lines, "")
        end
    end
    # Add whitespace and/or pipes to the left of the branches.
    pipe_centre = Int(floor((pipe_start + pipe_end) / 2))
    for idx in range(1, length=length(lines))
        pfx_start = "  "
        pfx_pipe = " "
        if idx >= pipe_start && idx <= pipe_end
            pfx_pipe = "|"
            if idx == pipe_centre
                # Introduce the current branch with /-, -- or \-.
                pfx_start = prefix
            end
        end
        lines[idx] = pfx_start * pfx_pipe * lines[idx]
    end

    return (lines, pipe_centre)
end

using DataFrames
using DataStructures

# A recursive function that takes PhyloTree as an argument
# Performs a single tree iteration. While walking from the top to bottom
# (since we don't have knowledge about leaves), saves the path. When reaching a
# leaf, iterates through already existing leaves (from the table) and calculates
# the distance between each pair, then writes results into the table.
# The distance is calculated by summing all paths until the common ancestor.
# The ancestor is extracted from using previously collected paths saved in the
# dictionary (by taking the intersection between two paths). Symmetric table
# entries are saved simultaneously. Returns the DataFrame table: df
function cophenetic(node::PhyloTree,path=PhyloTree[],dict=Dict{String,Array{PhyloTree}}(),df=DataFrame())
    #if the path is empty, we are at root, so start collecting path entries
    isempty(path) && push!(path,node)
    #if current node doesn't have children, the node is a leaf
    if isempty(node.children)
        #current path from root until the leaf is saved into the dictionary
        dict[node.name]=copy(path)
        column_length=length(names(df))
        values=Float64[]
        #for each already existing column (leaf)
        for name in names(df)
            #if the same node, the distance is 0
            if name == node.name
                push!(values,0.0)
            else
                #find the lowest common ancestor by getting paths of current pair from the dictionary
                #and take last element of their intersection
                common_ancestor=first(intersect(OrderedSet(reverse(dict[name])),OrderedSet(reverse(dict[node.name]))))
                cur=node
                sum1=0.0
                #sum of edges from the first path
                while cur != common_ancestor
                    sum1+=cur.dist
                    cur=cur.mother
                end
                cur=dict[name][end]
                sum2=0.0
                #sum of edges from the second path
                while cur != common_ancestor
                    sum2+=cur.dist
                    cur=cur.mother
                end
                push!(values,sum1+sum2)
            end
        end
        #fill table with symmetric distance values
        index=column_length+1
        push!(values,0.0)
        index!=1 && push!(df,values[begin:index-1])
        insertcols!(df,index,node.name=>values[begin:index])

    else
        #iterate though all children of current node
        for child in node.children
            #collecting path furter, recurively continue
            push!(path,child)
            cophenetic(child,path,dict,df)
            pop!(path)
        end
    end
    return df
end

#Function that finds the minimum value of the lower triangle in a symmetric matrix
function lower_triangle_min(matrix)
    min_value=Inf
    min_row=0
    min_col=0
    for i in 2:size(matrix)[1]
        for j in 1:i-1
            if matrix[i,j]<min_value
                min_value=matrix[i,j]
                min_row=i
                min_col=j
            end
        end
    end
    #returns a tuple(tuple(index of row, index of column), min value) 
    return (min_row,min_col),min_value
end


spanish = PhyloTree("Spanish", 1.7)
italian = PhyloTree("Italian", 1.7)
german = PhyloTree("German", 1.0)
dutch = PhyloTree("Dutch", 1.0)
english = PhyloTree("English", 1.5)
gothic = PhyloTree("Gothic", 0.5)
romance = PhyloTree("7", 2.3)
add_child!(romance, spanish)
add_child!(romance, italian)
cwg = PhyloTree("8", 0.5)
add_child!(cwg, dutch)
add_child!(cwg, german)
wg = PhyloTree("9", 1.0)
add_child!(wg, cwg)
add_child!(wg, english)
germanic = PhyloTree("10", 1.5)
add_child!(germanic, wg)
add_child!(germanic, gothic)
tree = PhyloTree("root", 0.0)
add_child!(tree, romance)
add_child!(tree, germanic)
print(tree)
show(cophenetic(tree))
