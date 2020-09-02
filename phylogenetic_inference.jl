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
