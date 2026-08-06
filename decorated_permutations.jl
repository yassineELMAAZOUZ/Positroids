"""
    decorated_permutations(n, k)

Return all decorated permutations of `1:n` having exactly `k` excedances.

A usual excedance is a position `i` for which `i < pi[i]`. Every fixed
point may have either a positive or a negative decoration, and a negatively
decorated fixed point contributes one excedance.

Encoding:
  * `i` at position `i` means a positively decorated fixed point;
  * `-i` at position `i` means a negatively decorated fixed point;
  * nonfixed entries are positive.

For example, `[-1, 3, 4, 2]` represents the permutation `[1, 3, 4, 2]`
with fixed point 1 negatively decorated. It has three excedances.
"""
function decorated_permutations(n::Integer, k::Integer)
    n >= 0 || throw(ArgumentError("n must be nonnegative"))
    0 <= k <= n || return Vector{Vector{Int}}()

    results = Vector{Vector{Int}}()
    decorated = Vector{Int}(undef, n)
    used = falses(n)

    function generate(position::Int, excedances::Int)
        excedances > k && return
        excedances + (n - position + 1) < k && return

        if position > n
            excedances == k && push!(results, copy(decorated))
            return
        end

        for value in 1:n
            used[value] && continue
            used[value] = true

            if value == position
                # Positively decorated fixed point: contributes 0.
                decorated[position] = value
                generate(position + 1, excedances)

                # Negatively decorated fixed point: contributes 1.
                decorated[position] = -value
                generate(position + 1, excedances + 1)
            else
                decorated[position] = value
                generate(position + 1, excedances + (position < value))
            end

            used[value] = false
        end
    end

    generate(1, 0)
    return results
end

"""Count excedances in an encoded decorated permutation."""
function decorated_excedances(p::AbstractVector{<:Integer})
    count(i -> i < abs(p[i]) || p[i] == -i, eachindex(p))
end

# Examples:
# decorated_permutations(4, 3)
# [-1, 3, 4, 2] in decorated_permutations(4, 3)  # true
# length(decorated_permutations(4, 3))
