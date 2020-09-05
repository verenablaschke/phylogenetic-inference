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

# Uses the neighbour joining algorithm to constructor a phylogenetic tree from
# a distance table.
# Note that even if the given cophenetic distance table allows for the
# construction of several different phylogenetic trees, this method returns
# only a single possible solution.
function nj(dist::DataFrame)
    dist = deepcopy(dist)
    # Keep track of all branches to assemble a full tree.
    branches = Dict{String, PhyloTree}(
        lang => PhyloTree(lang, -1.0) for lang in names(dist))

    while size(dist)[1] > 1
        langs = names(dist)
        n = length(langs)
        # Construct the lower triangle matrix Q.
        q = DataFrame(fill(Union{Float64, Missing}, n), langs, n)
        for i in 2:n
            for j in 1:i-1
                q[i, j] = (n - 2) * dist[i, j] - sum([dist[i, k] for k in range(1, stop=n)]) - sum([dist[j, k] for k in range(1, stop=n)])
            end
        end

        # The lowest value within Q shows the pair of languages that should be
        # merged next.
        merge_idx = lower_triangle_min(q)[1]
        merge_idx = sort([i for i in merge_idx])
        (f, g) = merge_idx

        # Determine the distance of the two branches that should be merged
        # to the node where they are joined together.
        dist_fg = dist[f, g]
        (lang1, lang2) = (langs[f], langs[g])
        (dist_f, dist_g) = (0.0, 0.0)
        for k in 1:length(langs)
            dist_f += dist[f, k]
            dist_g += dist[g, k]
        end
        if n > 2
            branches[lang1].dist = 0.5 * dist_fg + (1 / (2 * (n - 2))) * (dist_f - dist_g)
        else
            # Avoid dividing by 0 in the second part of the equation.
            # dist_f - dist_g = 0 in that case, anyway.
            branches[lang1].dist = 0.5 * dist_fg
        end
        branches[lang2].dist = dist_fg - branches[lang1].dist

        # Update the collection of tree branches.
        branch = PhyloTree(lang1 * lang2, 0.0)
        add_child!(branch, branches[lang1])
        add_child!(branch, branches[lang2])
        branches[lang1 * lang2] = branch
        delete!(branches, lang1)
        delete!(branches, lang2)

        # Find the distances from the newly joint node to all other branches.
        col = Array{Union{Missing, Float64},1}(undef, length(langs) + 1)
        col[1] = 0
        for k in 1:length(langs)
            if k in merge_idx
                continue
            end
            col[k + 1] = 0.5 * (dist[f, k] + dist[g, k] - dist_fg)
        end

        # Update the distance matrix: remove the now deprecated entries and add
        # the row and column vector of the joint node.
        deleteat!(col, [f + 1, g + 1])
        delete!(dist, merge_idx)  # remove rows
        select!(dist, Not(merge_idx))  # remove cols
        insertcols!(dist, 1, lang1 * lang2 => col[2:length(col)])
        foreach((c, v) -> insert!(c, 1, v), eachcol(dist), col)
    end

    # Only the complete phylogenetic tree is left in the branch collection
    # at this point.
    return collect(values(branches))[1]
end


# Returns a copy of the dataframe with (alphabetically) sorted column names
# and correspondingly updated cells.
function sort_dataframe(df::DataFrame)
    sorted_indices = sortperm(names(df))
    if sorted_indices == [i for i in 1:length(names(d))]
        return deepcopy(df)
    end
    sorted = df[!, sorted_indices]
    for i in 1:length(names(df))
        sorted[!, i] = sorted[!, i][sorted_indices]
    end
    return sorted
end


# Returns whether two dataframes are (approximately) identical.
function equivalent_tables(x::DataFrame, y::DataFrame)
    # isapprox returns a table with a 1 for each (approximately) identical cell
    # and a 0 for each different cell.
    for col in eachcol(isapprox.(sort_dataframe(x), sort_dataframe(y)))
        if 0 in col
            return false
        end
    end
    return true
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

println("The original tree:\n")
print(tree)
println("\n\nThe corresponding cophenetic table:\n")
distance_table = cophenetic(tree)
show(distance_table)
println("\n\nThe tree constructed via Neighbour Joining:\n")
nj_tree = nj(distance_table)
print(nj_tree)
println("\nThe cophenetic table inferred from this tree is identical to the original matrix:")
@show equivalent_tables(distance_table,  cophenetic(nj_tree))


# using RCall
# R"library(phangorn)"
# R"taxa <-c('English','Dutch','German','Italian','Spanish')"
# R"d <-as.dist(matrix(c(0.0,3.0,3.0,8.0,8.0,
#                        3.0,0.0,2.0,8.0,8.0,
#                        3.0,2.0,0.0,8.0,8.0,
#                        8.0,8.0,8.0,0.0,3.4,
#                        8.0,8.0,8.0,3.4,0.0),byrow=T,nrow=5,dimnames=list(taxa,taxa)))"
# R"nj.tree <- nj(d)"
# R"coph <- cophenetic(nj.tree)[taxa,taxa]"
# R"print(coph)"
# @rget coph
# @show coph
# @show coph == distance_table
