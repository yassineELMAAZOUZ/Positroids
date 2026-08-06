import Oscar
using LinearAlgebra
using Oscar: Perm, cycles, sign, young_tableau


function shiftedOrder(i::Int64,n::Int64)
    return vcat(i:n, 1:i-1)
end

# returns the cyclic interval a ..... b or a .... n 1 .... b-1 depending on a>b or b<a.
function cyclicInterval(n::Int64, i::Int64, j::Int64)
    if i<=j 
        return Vector(i:j)
    end

    if i>j 
        return reduce(vcat, (Vector(i:n), Vector(1:j-1)))
    end
end

# Compares two elements in the Gale order of i
function iCompare(i::Int64, n::Int64, x::Int64, y::Int64)
    
    L = shiftedOrder(i,n);
    
    posx = (findall(t->t==x, L))[1];
    posy = (findall(t->t==y, L))[1];

    return posx >= posy;
end


# Compares two Lists in the induced Gale order of i
function iCompareLists(i::Int64,n::Int64, X::Vector{Int64}, Y::Vector{Int64})
    
    L = shiftedOrder(i,n);
    
    Xind = [[(findall(t->t==x, L))[1], x] for x in X];
    Yind = [[(findall(t->t==y, L))[1], y] for y in Y];

    Xind = sort(Xind)
    Yind = sort(Yind)

    Xsor = [x[2] for x in Xind]
    Ysor = [y[2] for y in Yind]

    comparisons = unique([iCompare(i,n, Xsor[j], Ysor[j]) for j in 1:size(X)[1] ])

    if(size(comparisons)[1] > 1) 
        return false
    end
    
    if(comparisons[1] == false)
        return false
    end

    return true;
end


# Gets the minimal element in the ith Gale order of sequences.
function getMinList(i::Int64,n::Int64, B::Set{Vector{Int64}})

    X = collect(B)[1];

    for Y in B

        if iCompareLists(i,n, X, Y)
            X = Y;
        end

    end

    return X

end



# Gets the min Necklace given a positroid
function getNeckLace(n::Int64, B::Set{Vector{Int64}})

    return [getMinList(i,n,B) for i in 1:n]
    
end


function minGrassmannNecklace(n::Int64, B::Set{Vector{Int64}})
    return getNeckLace(n, B)
end


# Gets the minimal element in the ith Gale order of sequences.
function getMaxList(i::Int64,n::Int64, B::Set{Vector{Int64}})

    X = collect(B)[1];

    for Y in B

        if iCompareLists(i,n, Y, X)
            X = Y;
        end

    end

    return X

end


# Gets the max Necklace given a positroid
function getMaxNeckLace(n::Int64, B::Set{Vector{Int64}})

    return [getMaxList(i,n,B) for i in 1:n]
    
end


function maxGrassmannNecklace(n::Int64, B::Set{Vector{Int64}})
    return getMaxNeckLace(n, B)
end


function fromNecklaceToPositroid(k::Int64,n::Int64, N::Vector{Vector{Int64}})
    
    L = collect(Oscar.subsets(Vector(1:n), k));
    Intervals = [];
    Positroid = [];

    interval = [];

    for j in 1:n
        I = N[j];
        interval = [];

        for J in L
            
            if iCompareLists(j,n, J, I) 
                append!(interval, [J]);
            end
        end

        append!(Intervals, [interval]);

    end

    S = Set{Vector{Int64}}(Intervals[1]);
    
    for j in 2:n
        S = intersect(S, Set{Vector{Int64}}(Intervals[j]));
    end

    return S;
end



function fromNecklaceToDecoratedPermutation(N::Vector{Vector{Int64}})
    
    n = size(N)[1];
    p = Vector{Int64}(1:n);

    for i in 1:n
        
        j = i+1;

        if j == n+1
            j = 1;
        end

        I1 = Set(N[i])
        I2 = Set(N[j])

        if I1 != I2
            p_i = collect(setdiff(I2,I1))[1]
            p[p_i] = i
        end
        
        if I1 == I2

            if i in I1
                p[i] = -i
            end
            if ! (i in I1)
                p[i] = i
            end

        end
    end
    return p
end


# Computes the Grassmann Necklace associated to the decorated permutation
function fromDecoratedPermToNecklace(p::Vector{Int64})

    n = size(p)[1];

    N = Vector{Vector{Int64}}([]);
    
    for k in 1:n

        I = Vector{Int64}([]);

        for i in 1:n
            
            if abs(p[i]) == i
                if abs(p[i]) == - p[i]
                    append!(I, [i])
                end
            else 
                if iCompare(k,n,   abs(p[i]) , i)
                    append!(I, [i])
                end
            end
        end

        append!(N, [I]);
    end


    return N;
    
end


function fromDecoratedPermToMinNecklace(p::Vector{Int64})
    return fromDecoratedPermToNecklace(p)
end


function minGrassmannNecklace(p::Vector{Int64})
    return fromDecoratedPermToMinNecklace(p)
end


function fromDecoratedPermToMaxNecklace(p::Vector{Int64})

    _validate_decorated_permutation(p)

    n = length(p)
    k = countExceedences(p)

    if k == 0
        return [Int64[] for i in 1:n]
    end

    if k == n
        return [Vector(1:n) for i in 1:n]
    end

    return getMaxNeckLace(n, fromDecoratedPermToPositroid(k, n, p))
end


function maxGrassmannNecklace(p::Vector{Int64})
    return fromDecoratedPermToMaxNecklace(p)
end


function _affine_value_from_decorated_perm(p::Vector{Int64}, i::Int64)

    n = length(p)
    pi = abs(p[i])

    if pi == i
        return p[i] < 0 ? i + n : i
    end

    return pi > i ? pi : pi + n
end


function reverseGrassmannNecklace(p::Vector{Int64})

    _validate_decorated_permutation(p)

    n = length(p)
    shifts = [_affine_value_from_decorated_perm(p, i) - i for i in 1:n]
    N = Vector{Vector{Int64}}([])

    for q in 1:n
        I = Vector{Int64}([])

        for r in 1:n
            representative = r <= q ? r : r - n

            if representative <= q && q < representative + shifts[r]
                push!(I, modN(representative, n))
            end
        end

        push!(N, sort(I))
    end

    return N
end


# This is either wrong or the next function is wrong
function fromDecoratedPermToPositroid(k::Int64,n::Int64, p::Vector{Int64})
    
    n = size(p)[1];
    N = fromDecoratedPermToNecklace(p);

    Pos = fromNecklaceToPositroid(k,n, N); 
    return Pos;
end


function fromPositroidToPermutation(n::Int64, M::Set{Vector{Int64}})
    
    N = getNeckLace(n,M);

    return fromNecklaceToDecoratedPermutation(N);
end


# Function to generate vertices of a regular n-gon
function generate_ngon_vertices(n::Int, r::Float64 = 1.0)
    angles = range(0, 2π, length=n+1)[1:end-1]  # n angles evenly spaced
    return [(r*cos(-θ-π/2), r*sin(-θ-π/2)) for θ in angles]
end


# Function to draw the n-gon and chords
function draw_chords(p::Vector{Int64})

    _load_plotting()

    r = 1.0;

    n = size(p)[1]

    vertices = generate_ngon_vertices(n)

    x_coords = [v[1] for v in vertices]
    y_coords = [v[2] for v in vertices]


    # Plot the disk (circle)
    θ = range(0, 2π, length=100)
    Plots.plot(r*cos.(θ), r*sin.(θ), seriestype = :path, lw=2, label="Disk", legend=false, ratio=1, fill=:blue, fillalpha=0.1)

    # Plot the n-gon
    #plot!(x_coords, y_coords, seriestype = :shape, lw=3, label="n-gon", legend=false, ratio=1)

    # Draw all possible chords
    for i in 1:n
        
        
        if abs(p[i])!= i
            GR.setarrowsize(2)
            j = abs(p[i]);
            Plots.plot!([x_coords[i], x_coords[j]], [y_coords[i], y_coords[j]], lw=4, linecolor=:black, label=false, arrow=(:closed, 2.0))
        end


        if  p[i] < 0
            GR.setarrowsize(2)
            x_text = x_coords[i] * 1.2   # Slightly offset the text position
            y_text = y_coords[i] * 1.2
            Plots.annotate!(x_text, y_text, Plots.text("-", "red", 26))
        end

    end

    # Number the vertices clockwise from 1 to n
    for i in 1:n
       
        label = string(i)
        x_text = x_coords[i] * 1.05  # Slightly offset the text position
        y_text = y_coords[i] * 1.05
        Plots.annotate!(x_text, y_text, label)
       
       
        #if i%2 == 1
        #    x_text = x_coords[i] * 1.2  # Slightly offset the text position
        #    y_text = y_coords[i] * 1.2
        #    annotate!(x_text, y_text, text("+", "blue", 20))
        #end
        
        #if i%2 == 0
        #    x_text = x_coords[i] * 1.2  # Slightly offset the text position
        #    y_text = y_coords[i] * 1.2
        #    annotate!(x_text, y_text, text("-", "red", 20))
        #end
    end


    # Display the plot
    return Plots.plot!(showaxis = false, ticks = false, ylimits=(-1.5,1.5),  xlimits=(-1.5,1.5))
end


















function signOfList(I::Vector{Int64})

    indexedI = [[I[i], i] for i in 1:size(I)[1]]

    sortedI = sort(indexedI);

    p = [x[2] for x in sortedI];

    return sign(Perm(p));
end

function is_orthogonal(k::Int64, n::Int64, M::Set{Vector{Int64}}, verbose=false)

    W = collect(Oscar.subsets(Vector(1:n), k-1));

    N = binomial(n,k-1);

    for i in 1:N
        
        for j in i:N

            I = W[i];
            J = W[j];

            IJcomp = sort(collect(setdiff( Set(Vector(1:n)) ,union(Set(I), Set(J)) )));

            L1 = filter( l-> signOfList(reduce(vcat, (I, [l])) ) * signOfList(reduce(vcat, (J, [l])) ) * ( sort(reduce(vcat, (I, [l])) ) in M) * ( sort(reduce(vcat, (J, [l])) ) in M) * (-1)^l ==  1 , IJcomp) 
            L2 = filter( l-> signOfList(reduce(vcat, (I, [l])) ) * signOfList(reduce(vcat, (J, [l])) ) * ( sort(reduce(vcat, (I, [l])) ) in M) * ( sort(reduce(vcat, (J, [l])) ) in M) * (-1)^l == -1 , IJcomp) 

            if  ( isempty(L1) &&  ! isempty(L2) ) || (!isempty(L1) && isempty(L2)) 
                
                if verbose
                    println("Problem at: I= ", I," and J= ",J);
                    println("+ side indices : ",L1, " and - side indices:", L2 )
                end

                return false;

            end

        end

    end


    return true
end


function is_orthogonal(p::Vector{Int64}, verbose=false)
    k = countExceedences(p);
    n = size(p)[1];
    M = fromDecoratedPermToPositroid(k,n,p)
    return is_orthogonal(k,n,M);
end

function alignmentNumber(p::Vector{Int64})
    
    n = size(p)[1]

    counter = 0


    for i in 1:n
        for j in 1:n

            if i != j


                pi = abs(p[i])
                pj = abs(p[j])

                

                i_pj_interval = cyclicInterval(n, i,pj)
                pj_i_interval = cyclicInterval(n, pj,i)
                
                v = (pi in i_pj_interval) &&  (j in pj_i_interval)

                if (pi == i)
                    v = v && (p[i] < 0)
                end

                if (pj == j)
                    v = v && (p[j] > 0)
                end
                
                if v 
                    
                    counter = counter + 1;
                end

            end

        end
    end
    
    return counter
    
end


function dimensionOfPositroid(k::Int64, n::Int64, M::Set{Vector{Int64}})

    p = fromPositroidToPermutation(n, M);

    return k*(n-k) - alignmentNumber(p);
    
end

function dimensionOfPermutation(k::Int64, n::Int64, p::Vector{Int64})
    
    return k*(n-k) - alignmentNumber(p);
    
end


function numVertices(n)

    if n in 0:3
        return 0
    end

    s = 0

    for l in 1:trunc(Int64, n/2)

        if n%2 ==0

            s = s+ (l-1)^2 +   trunc(Int64, ((n-2*l)/2)^2);

        end

        if n%2 == 1

            s = s+ (l-1)^2 +   trunc(Int64, ((n+1-2*l)/2) * ((n-1-2*l)/2) );

        end

    end
    
    return s + numVertices(n-1)
    
end

function countExceedences(p::Vector{Int64})
    return sum([( i < abs(p[i])) + (i == - p[i]) for i in 1:size(p)[1]]);
end



function permutationFromCycles(n::Int64, C::Vector{Vector{Int64}})

    p = [0 for i in 1:n]

    for i in 1:n
        
        for I in C
            if i in I
                j = findall(t->t==i, I)[1];

                if j < size(I)[1]
                    p[i] = I[j+1]
                else
                    p[i] = I[1]
                end
            end
        end


        if p[i] == 0
            p[i] = i
        end
    end

    return p;


end




# Define a function to check the orientation of three points
function orientation(p, q, r)
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
    if val > 0
        return 1  # Clockwise
    elseif val < 0
        return -1 # Counter-clockwise
    else
        return 0  # Collinear
    end
end

# Check if a point is on the segment defined by two other points
function on_segment(p, q, r)
    return min(p[1], r[1]) <= q[1] <= max(p[1], r[1]) && min(p[2], r[2]) <= q[2] <= max(p[2], r[2])
end

# Check if two line segments (p1, q1) and (p2, q2) intersect
function do_intersect(p1, q1, p2, q2)
    # Find the orientations of the ordered triplets (p1, q1, p2), (p1, q1, q2), (p2, q2, p1), (p2, q2, q1)
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case: segments intersect if orientations are different
    if o1 != o2 && o3 != o4
        return true
    end

    # Special cases
    if o1 == 0 && on_segment(p1, p2, q1)
        return true
    end
    if o2 == 0 && on_segment(p1, q2, q1)
        return true
    end
    if o3 == 0 && on_segment(p2, p1, q2)
        return true
    end
    if o4 == 0 && on_segment(p2, q1, q2)
        return true
    end

    return false
end

# Function to check if two polygons intersect
function polygons_intersect(poly1, poly2)
    n1 = length(poly1)
    n2 = length(poly2)

    # Check every edge of poly1 against every edge of poly2
    for i in 1:n1
        p1 = poly1[i]
        q1 = poly1[(i % n1) + 1]
        for j in 1:n2
            p2 = poly2[j]
            q2 = poly2[(j % n2) + 1]
            if do_intersect(p1, q1, p2, q2)
                return true  # An intersection was found
            end
        end
    end
    return false  # No intersection found
end

# Function to group only the data based on polygon connected components
function group_data_by_polygon_components(polygons_with_data)
    n = length(polygons_with_data)
    visited = falses(n)  # Array to keep track of visited polygons
    components = []      # Array to store connected components of data

    # Define a DFS function to explore connected polygons and collect data
    function dfs(polygon_idx, component_data)
        visited[polygon_idx] = true
        # Collect the data associated with the polygon
        push!(component_data, polygons_with_data[polygon_idx][2])

        # Check all other polygons for intersection
        for i in 1:n
            if !visited[i] && polygons_intersect(polygons_with_data[polygon_idx][1], polygons_with_data[i][1])
                dfs(i, component_data)
            end
        end
    end

    # Loop over all polygons and perform DFS for unvisited polygons
    for i in 1:n
        if !visited[i]
            component_data = []
            dfs(i, component_data)
            push!(components, component_data)
        end
    end

    return components
end


function orthogonal_dimension(k::Int64, n::Int64, M::Set{Vector{Int64}})

 

    p = fromPositroidToPermutation(n,M);
    k = countExceedences(p);
    vertices = generate_ngon_vertices(n);


    perm_p = Perm(p);
    cyc = filter( c->size(c)[1] > 1 ,collect(cycles(perm_p)));
    polygons_with_cycles = [ [[vertices[i] for i in c], c] for c in cyc];




    N = 0;

    for i in 1:size(cyc)[1]
    
        ki = countExceedences(permutationFromCycles(n,[polygons_with_cycles[i][2]]));

        N = N + binomial(ki+1,2);

        for j in i+1:size(cyc)[1]
            if polygons_intersect(polygons_with_cycles[i][1], polygons_with_cycles[j][1])
                
                kj = countExceedences(permutationFromCycles(n,[polygons_with_cycles[j][2]]))
                N = N + ki*kj
            end
        end
    end


    return dimensionOfPositroid(k,n,M) - N;

end





function orthogonal_dimensionOfPermutation(p::Vector{Int64})
    n = size(p)[1];
    k = countExceedences(p);

    M = fromDecoratedPermToPositroid(k,n,p)
    
    return orthogonal_dimension(k,n,M);
    
end







function inversionsOfPerm(p::Vector{Int64})

    n = size(p)[1];

    Invs = []

    for i in 1:n-1
        for j in i+1:n

            if p[i]> p[j]
                append!(Invs, [(i,j)]);
            end
        end
    end
    
    return Invs
end


function pminversionsOfPerm(p::Vector{Int64})

    n = size(p)[1];

    Invs = []

    for i in 1:n-1
        for j in i+1:n

            if p[i]> p[j] && (-1)^(i+j) == -1
                append!(Invs, [(i,j)]);
            end
        end
    end
    
    return size(Invs)[1]
end



# The order on the positroids is not containement!
function is_codim1_face(p::Vector{Int64}, q::Vector{Int64})
    
    n = size(p)[1]
    k = countExceedences(p);
    
    if countExceedences(q) !=k
        return "The two permutations should have the same k";
    end

    d = orthogonal_dimensionOfPermutation(p)

    if d != orthogonal_dimensionOfPermutation(q) -1
        return false;
    end

    Mp = fromDecoratedPermToPositroid(k,n,p)
    Mq = fromDecoratedPermToPositroid(k,n,q)

    return isempty(setdiff(Mp, Mq));
end


function is_child(p::Vector{Int64}, q::Vector{Int64})
    
    n = size(p)[1]
    k = countExceedences(p);
    
    if countExceedences(q) !=k
        return "The two permutations should have the same k";
    end

    Mp = fromDecoratedPermToPositroid(k,n,p)
    Mq = fromDecoratedPermToPositroid(k,n,q)

    return isempty(setdiff(Mp, Mq));
end

function _validate_decorated_permutation(p::Vector{Int64})

    n = length(p)

    if any(x -> x == 0 || abs(x) > n, p)
        error("A decorated permutation of [1,$n] can only contain signed values in 1:$n.")
    end

    if sort(abs.(p)) != Vector(1:n)
        error("The absolute values of a decorated permutation must be a permutation of 1:$n.")
    end

    for i in 1:n
        if abs(p[i]) != i && p[i] < 0
            error("Only fixed points may carry a negative decoration; p[$i] = $(p[i]).")
        end
    end

    return true
end


function _le_shape_data(I::Vector{Int64}, n::Int64)

    I = sort(I)
    k = length(I)
    Ic = sort(collect(setdiff(Set(Vector(1:n)), Set(I))), rev=true)

    row_lengths = [n - k + r - I[r] for r in 1:k]

    if any(l -> l < 0 || l > n-k, row_lengths)
        error("The first necklace term $I does not determine a Young diagram inside a $k x $(n-k) rectangle.")
    end

    if !issorted(row_lengths, rev=true)
        error("The row lengths $row_lengths are not a partition shape.")
    end

    return I, Ic, row_lengths
end


function le_diagram_matrix(p::Vector{Int64})

    _validate_decorated_permutation(p)

    n = length(p)
    N = fromDecoratedPermToNecklace(p)
    I, Ic, row_lengths = _le_shape_data(N[1], n)

    T = [zeros(Int, row_lengths[r]) for r in 1:length(row_lengths)]

    row_index = Dict(I[r] => r for r in 1:length(I))
    column_index = Dict(Ic[c] => c for c in 1:length(Ic))

    for j in 2:n

        J = sort(N[j])

        A = sort(collect(setdiff(Set(I), Set(J))), rev=true)
        B = sort(collect(setdiff(Set(J), Set(I))))

        if length(A) != length(B)
            error("Invalid Grassmann necklace step at index $j: removed $A but added $B.")
        end

        for l in 1:length(A)
            a, b = A[l], B[l]

            rrow = row_index[a]
            ccolumn = column_index[b]

            if ccolumn > length(T[rrow])
                error("The box labeled by row $a and column $b lies outside the Young diagram.")
            end

            T[rrow][ccolumn] = 1
        end
    end

    return T
end


function le_diagram(p::Vector{Int64}; as_tableau::Bool=true)

    T = le_diagram_matrix(p)

    if !as_tableau || !isdefined(@__MODULE__, :young_tableau)
        return T
    end

    try
        return young_tableau(T)
    catch err
        if err isa ArgumentError && all(isempty, T)
            return T
        end

        rethrow()
    end
end



function is_Omega_symplectic(k::Int64, n::Int64, M::Set{Vector{Int64}})

    m = n//2;

    V = collect(Oscar.subsets(n, k-2));

    for S in V 

        A_odd  = filter(t-> ((t + size(  filter( i-> i>t, S))[1]+size(  filter( i-> i>t+m, S))[1] ) % 2 == 1   && sort(unique(vcat(S,[t,m+t]))) in M)  ,   Vector(1:m) );
        A_even = filter(t-> ((t + size(  filter( i-> i>t, S))[1]+size(  filter( i-> i>t+m, S))[1] ) % 2 == 0   && sort(unique(vcat(S,[t,m+t]))) in M)  ,   Vector(1:m) );

        if  ( isempty(A_odd) &&  ! isempty(A_even) ) || (!isempty(A_odd) && isempty(A_even)) 
          return false
        end

    end

    return true
end




function is_Omega_symplecticTwo(k::Int64, n::Int64, M::Set{Vector{Int64}})

    m = n//2;

    V = collect(Oscar.subsets(n, k-2));

    for S in V 

        A_odd  = filter(t-> ((t + size(  filter( i-> i>t, S))[1]+size(  filter( i-> i>t+m, S))[1] ) % 2 == 1   && sort(unique(vcat(S,[t,m+t]))) in M)  ,   Vector(1:m) );
        A_even = filter(t-> ((t + size(  filter( i-> i>t, S))[1]+size(  filter( i-> i>t+m, S))[1] ) % 2 == 0   && sort(unique(vcat(S,[t,m+t]))) in M)  ,   Vector(1:m) );


        if (size(A_odd)[1] in [0,1]) && (size(A_even)[1] in [0,1])

            if  (isempty(A_odd) &&  ! isempty(A_even)) || (!isempty(A_odd) && isempty(A_even)) 
                return false
            end

        end

    end

    return true
end




function symplectic_Omega_dimension(k::Int64, n::Int64, M::Set{Vector{Int64}})

    if k > 2 return "Not implemented for k>2 yet!" end 

    m = n//2;

    A_odd  = filter(t-> (t%2==1 && [t,m+t] in M) == true, Vector(1:m));
    A_even = filter(t-> (t%2==0 && [t,m+t] in M) == true, Vector(1:m));

    if !isempty(A_odd) && !isempty(A_even)
        return dimensionOfPositroid(k,n,M)-1;
    end

    return dimensionOfPositroid(k,n,M);
end




function modN(i,N)

    if i > N
        return i - N
    elseif i < 1
        return i + N
    else
        return i
    end
    
end




function is_E_symplectic(k::Int64, n::Int64, M::Set{Vector{Int64}})

    m = n//2;

    V = collect(Oscar.subsets(n, k-2));

    for S in V 

        A_odd  = filter(t-> (t%2==1 && sort(unique(vcat(S,[t,2*m-t+1]))) in M) == true, Vector(1:m));
        A_even = filter(t-> (t%2==0 && sort(unique(vcat(S,[t,2*m-t+1]))) in M) == true, Vector(1:m));

        if  ( isempty(A_odd) &&  ! isempty(A_even) ) || (!isempty(A_odd) && isempty(A_even)) 
            return false
        end

    end

    return true
end



function symplectic_E_dimension(k::Int64, n::Int64, M::Set{Vector{Int64}})

    if k > 2 return "Not implemented for k>2 yet!" end 

    m = n//2;

    A_odd  = filter(t-> (t%2==1 && [t,2*m-t+1] in M) == true, Vector(1:m));
    A_even = filter(t-> (t%2==0 && [t,2*m-t+1] in M) == true, Vector(1:m));

    if !isempty(A_odd) && !isempty(A_even)
        return dimensionOfPositroid(k,n,M)-1;
    end

    return dimensionOfPositroid(k,n,M);
end



function up(k::Int64, n::Int64,  M::Set{Vector{Int64}})

    if n!=2k+1 
        error("n != 2k + 1");
    end

    newM = Set{Vector{Int64}}([]);

    for I in M 
        J = reduce(vcat, [I, [n+1]])
        push!(newM, J);

        Jc = filter( t-> ! (t in J) == true,  1:n+1);
        push!(newM, Jc);
    end

    return newM;
end


function randomMatching(n::Int64)

    if n%2 != 0
        error("n is not even!!");
    end

    V = Vector(1:n);

    C = Vector{Vector{Int64}}([]);
    
    for i in 1:n//2
    
        a = rand(V);
        V = filter(j-> !(j==a)==true, V);
        b = rand(V);
        V = filter(j-> !(j==b)==true, V);
    
        C = append!(C, [[a,b]]);
    end
    

    return permutationFromCycles(n, C);

end












function twistPerm(p)

    m = size(p)[1];
    n = trunc(Int64, m/2)
    q = Vector{Int64}([]);

    pp = abs.(p);

    for i in 1:m

        e = 1;

        if pp[i] == -p[i]
            e = -1;
        end

        push!(q, e* modN( pp[i]+n,m));

    end


    return q

end





function untwistPerm(p)

    m = size(p)[1];
    n = trunc(Int64, m/2)
    q = Vector{Int64}([]);

    pp = abs.(p);


    for i in 1:m
        
        e = 1

        if p[i] < 0
            e = -1;
        end

        push!(q,  e * modN( pp[i] - n, m));
    end

    return q

end

function twistNecklace(I:: Vector{Vector{Int64}})

    m = size(I)[1];
    n = trunc(Int64,m/2);
    
    Res = Vector{Vector{Int64}}([]);

    J = [I[ modN(i+n,m) ] for i in 1:m];

    for i in 1:m

        S = J[i];

        rotS = [modN(x+n,m) for x in S]

        push!(Res, filter( x-> !(x in rotS) , 1:m))

    end

    return Res

end


function _validate_twist_necklace(N::Vector{Vector{Int64}}, k::Int64, n::Int64)

    if length(N) != n
        error("Expected a Grassmann necklace with $n terms, but got $(length(N)).")
    end

    for i in 1:n
        I = sort(N[i])

        if length(I) != k
            error("Necklace term $i has length $(length(I)); expected $k.")
        end

        if length(unique(I)) != k || any(j -> j < 1 || j > n, I)
            error("Necklace term $i is not a $k-subset of 1:$n: $(N[i]).")
        end
    end

    return [sort(I) for I in N]
end


function _is_zero_for_twist(x, atol)
    if x isa AbstractFloat
        return abs(x) <= atol
    end

    return iszero(x)
end


function _is_nonzero_for_twist(x, atol)
    if x isa AbstractFloat
        return abs(x) > atol
    end

    return !iszero(x)
end


function _check_twist_basis(A::AbstractMatrix, I::Vector{Int64}, i::Int64, atol)

    d = det(A[:, I])

    if !_is_nonzero_for_twist(d, atol)
        error("The columns indexed by necklace term I_$i = $I are not a basis for the given matrix.")
    end

    return nothing
end


function _check_twist_loop_column(A::AbstractMatrix, i::Int64, atol)

    if any(x -> !_is_zero_for_twist(x, atol), A[:, i])
        error("Necklace term I_$i does not contain $i, so the cell expects column $i to be a loop/zero column.")
    end

    return nothing
end


function _solve_twist_column(A::AbstractMatrix, I::Vector{Int64}, position::Int64)

    rhs = zeros(eltype(A), size(A, 1))
    rhs[position] = one(eltype(A))

    return transpose(A[:, I]) \ rhs
end


# Right twist of a point A in the positroid cell with Grassmann necklace N.
# This is symbolic-friendly: pass a matrix over SymPy.Sym, rationals, or another
# exact field/ring and the returned entries are computed by exact linear solves.
function right_twist(A::AbstractMatrix, N::Vector{Vector{Int64}}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)
    N = _validate_twist_necklace(N, k, n)

    twisted_columns = Vector{Any}(undef, n)

    for i in 1:n
        I = N[i]
        position = findfirst(==(i), I)

        if isnothing(position)
            if check_cell
                _check_twist_loop_column(A, i, atol)
            end

            twisted_columns[i] = nothing
            continue
        end

        if check_cell
            _check_twist_basis(A, I, i, atol)
        end

        twisted_columns[i] = _solve_twist_column(A, I, position)
    end

    nonloop_column = findfirst(c -> !isnothing(c), twisted_columns)

    if isnothing(nonloop_column)
        return zeros(eltype(A), k, n)
    end

    zero_column = zero(twisted_columns[nonloop_column])

    for i in 1:n
        if isnothing(twisted_columns[i])
            twisted_columns[i] = zero_column
        end
    end

    return hcat(twisted_columns...)
end


function right_twist(A::AbstractMatrix, p::Vector{Int64}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)

    if length(p) != n
        error("The decorated permutation has length $(length(p)), but the matrix has $n columns.")
    end

    if countExceedences(p) != k
        error("The decorated permutation has k = $(countExceedences(p)), but the matrix has $k rows.")
    end

    return right_twist(A, fromDecoratedPermToNecklace(p); check_cell=check_cell, atol=atol)
end


function right_twist(A::AbstractMatrix, M::Set{Vector{Int64}}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)

    if !isempty(M) && length(collect(M)[1]) != k
        error("The positroid consists of $(length(collect(M)[1]))-subsets, but the matrix has $k rows.")
    end

    return right_twist(A, getNeckLace(n, M); check_cell=check_cell, atol=atol)
end


function left_twist(A::AbstractMatrix, N::Vector{Vector{Int64}}; check_cell::Bool=true, atol=1e-10)

    return right_twist(A, N; check_cell=check_cell, atol=atol)
end


function left_twist(A::AbstractMatrix, p::Vector{Int64}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)

    if length(p) != n
        error("The decorated permutation has length $(length(p)), but the matrix has $n columns.")
    end

    if countExceedences(p) != k
        error("The decorated permutation has k = $(countExceedences(p)), but the matrix has $k rows.")
    end

    return left_twist(A, reverseGrassmannNecklace(p); check_cell=check_cell, atol=atol)
end


function left_twist(A::AbstractMatrix, M::Set{Vector{Int64}}; check_cell::Bool=true, atol=1e-10)

    k, n = size(A)

    if !isempty(M) && length(collect(M)[1]) != k
        error("The positroid consists of $(length(collect(M)[1]))-subsets, but the matrix has $k rows.")
    end

    return left_twist(A, reverseGrassmannNecklace(fromPositroidToPermutation(n, M)); check_cell=check_cell, atol=atol)
end


function twist_map(A::AbstractMatrix, cell; side::Symbol=:right, check_cell::Bool=true, atol=1e-10)

    if side == :right
        return right_twist(A, cell; check_cell=check_cell, atol=atol)
    end

    if side == :left
        return left_twist(A, cell; check_cell=check_cell, atol=atol)
    end

    error("Unknown twist side $side. Use side=:right or side=:left.")
end





## Function that computes dimension of a decorated involution


function dimensionOfLGperm(p::Vector{Int64})

    q = abs.(twistPerm(p));

    m = size(p)[1];
    n = trunc(Int64, m/2);

    chords = Vector{Vector{Int64}}([]);
    
    visited = Vector{Int64}([]);

    for i in 1:m

        if !(i in visited)

            j = q[i];

            if j != i

                push!(chords, [i,j])
                append!(visited, [i,j]);   
            end

            append!(visited, [i,j]);   
        end

    end


    D = 0;

    for c in chords
        
        D = D + min( modN(c[1] - c[2],m), modN(c[2] - c[1],m)  );

    end

    if size(chords)[1] > 1   

            for a in 1:size(chords)[1]-1

                for b in a+1:size(chords)[1]

                    c1 = chords[a];
                    c2 = chords[b];

                    i,j = 1, modN(c1[2] - c1[1]+1,m)
                    k,l = modN(c2[1] - c1[1]+1,m) , modN(c2[2] - c1[1]+1,m)    
                    
                    if min(k,l) < j && j < max(k,l) 

                        D = D-1;

                    end

                end
            end
    end


    return trunc(Int64, n*(n+1)/2) - D

end










function iota(I::Vector{Int64})

    n = size(I)[1];
    J = [modN(i+n,2*n) for i in I];    
    return sort(filter( x-> !(x in J) , 1:2*n))

end





function eta(i::Int64,n::Int64)

    modN(i+1,2*n)
end


function eta(I::Vector{Int64},n::Int64)
    return sort([eta(i,n) for i in I])
end




function iota2(I::Vector{Int64})

    n = size(I)[1];
    J = eta(I,n)
    return sort(filter( x-> !(x in J) , 1:2*n))

end





function is_cyclo_isotropic(k::Int64, n::Int64, M::Set{Vector{Int64}})

    for I in M

        if !(iota2(I) in M)
            return false
        end
    end
    return true
end
