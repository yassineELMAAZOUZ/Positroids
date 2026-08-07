using Test
using Positroids
include("screenshot_reduction_fixture.jl")

@testset "decorated permutations" begin
    @test decorated_permutations(0, 0) == [Int[]]
    @test decorated_permutations(3, -1) == Vector{Vector{Int}}()
    @test [-1, 3, 4, 2] in decorated_permutations(4, 3)
    @test decorated_excedances([-1, 3, 4, 2]) == 3
    @test countExceedences([-1, 3, 4, 2]) == 3
    @test decorated_excedances([-1,2,3,4])==1 # -1 is a coloop
    @test decorated_excedances([1,2,3,4])==0  # positive fixed points are loops
    @test Positroids._lollipop_color(-1,1)==:white
    @test Positroids._lollipop_color(1,1)==:black
    @test Positroids._lollipop_decoration(:white,1)==-1
    @test Positroids._lollipop_decoration(:black,1)==1
    @test positroid_rank([4,3,1,5,2])==2
    @test positroid_rank([3,5,2,1,4];permutation_convention=:source)==2
    @test_throws ArgumentError positroid_rank([4,3,1,5,2];
                                               permutation_convention=:unknown)
end

@testset "necklace/permutation round trips" begin
    for n in 1:5, k in 0:n, p in decorated_permutations(n, k)
        N = fromDecoratedPermToNecklace(p)
        @test all(length(I) == k for I in N)
        @test fromNecklaceToDecoratedPermutation(N) == p
        M = fromDecoratedPermToPositroid(k, n, p)
        @test fromPositroidToPermutation(n, M) == p
    end
    @test fromDecoratedPermToPositroid(2,4,[-1,-2,3,4])==Set([[1,2]])
    @test fromDecoratedPermToPositroid(0,4,[1,2,3,4])==Set([Int[]])
end

@testset "validation" begin
    @test_throws ArgumentError fromDecoratedPermToPositroid(1, 4, [2, 1])
    @test_throws ArgumentError fromDecoratedPermToPositroid(0, 2, [2, 1])
end

@testset "Le diagrams and twists" begin
    p = [3, 4, 1, 2]
    @test le_diagram_matrix(p) isa Vector
    A = [1.0 0.0 1.0 1.0; 0.0 1.0 1.0 2.0]
    T = right_twist(A, p; check_cell=false)
    @test size(T) == size(A)
end

@testset "Plucker positivity" begin
    A = [1 0 -1; 0 1 1]
    @test plucker_coordinates(A) == Dict((1, 2) => 1, (1, 3) => 1, (2, 3) => 1)
    @test is_totally_positive(A)
    @test is_totally_nonnegative([1 0 -1; 0 1 0])
    @test !is_totally_positive([1 0 -1; 0 1 0])
    @test !is_totally_nonnegative([1 0 1; 0 1 1])
end

@testset "plabic bridge graphs" begin
    for n in 1:5, k in 0:n, p in decorated_permutations(n, k)
        G = plabic_graph(p)
        @test G.n == n
        @test trip_permutation(G) == p
        degrees = Positroids._plabic_degrees(G)
        @test all(degrees[1:n] .== 1)
        @test all(d > 0 && d != 2 for d in degrees[(n + 1):end])
        @test all(c in (:boundary, :black, :white) for c in G.colors)
        for i in 1:n
            abs(p[i])==i || continue
            neighbor=only([u==i ? v : u for (u,v) in G.edges if u==i || v==i])
            @test G.colors[neighbor]==(p[i]<0 ? :white : :black)
        end
        @test is_bipartite(G)
        @test is_reduced(G)
        @test validate_plabic_graph(G)
    end
end

@testset "interactive graph transforms" begin
    p=[4,3,1,5,2]
    G=plabic_graph(p)
    shift=2
    H=Positroids._cyclically_relabel_plabic_graph(G,shift)
    relabel(v)=v<=G.n ? mod1(v+shift,G.n) : v
    expected_edges=sort(unique([(min(relabel(a),relabel(b)),max(relabel(a),relabel(b)))
                                for (a,b) in G.edges]))
    @test H.edges==expected_edges
    @test all(H.positions[mod1(i+shift,G.n)]==G.positions[i] for i in 1:G.n)
    @test is_bipartite(H)
    @test is_reduced(H)
    @test last.(graph_trips(plabic_embedding(H;iterations=500,restarts=4)))==
          abs.(H.permutation)
    restored=Positroids._cyclically_relabel_plabic_graph(H,-shift)
    @test restored.permutation==G.permutation
    @test restored.colors==G.colors
    @test restored.edges==sort(G.edges)
    @test restored.positions==G.positions

    swapped=Positroids._swap_plabic_colors(G)
    @test swapped.permutation==[3,5,2,1,4]
    @test swapped.colors==[color==:black ? :white : color==:white ? :black : color
                            for color in G.colors]
    @test is_bipartite(swapped)
    @test is_reduced(swapped)
    @test last.(graph_trips(plabic_embedding(swapped;iterations=500,restarts=4)))==
          abs.(swapped.permutation)
    swapped_back=Positroids._swap_plabic_colors(swapped)
    @test swapped_back.permutation==G.permutation
    @test swapped_back.colors==G.colors
    @test swapped_back.edges==G.edges
    @test swapped_back.positions==G.positions
    @test Positroids._swap_plabic_colors(plabic_graph([-1,-2,3,4])).permutation==
          [1,2,-3,-4]
end


@testset "bivalent reduction" begin
    @test PACKAGE_BUILD == v"1.3.0"
    G = plabic_graph([2, 1])
    @test length(G.colors) == 2
    @test G.edges == [(1, 2)]
    @test is_reduced(G)
    @test vertex_degrees(G) == [1, 1]
    @test G.positions[1][1] ≈ 0.0 atol=1e-12
    @test G.positions[1][2] ≈ -1.0 atol=1e-12
    rotation=Positroids._rotation_from_positions(G,unique(vcat(G.edges,[(1,2),(2,1)])))
    @test rotation[1] == [2,2,2]
    @test rotation[2] == [1,1,1]
    rebuilt_direct=plabic_graph_from_drawing(2,Symbol[],[(1,2)],Tuple{Float64,Float64}[])
    @test trip_permutation(rebuilt_direct)==[2,1]
    @test is_reduced(rebuilt_direct)
end

@testset "graphs from disk drawings" begin
    colors=[:white,:white,:black,:black]
    positions=[(0.0,-0.8),(-0.8,0.0),(0.0,0.8),(0.8,0.0)]
    edges=[(1,5),(2,6),(3,7),(4,8)]
    G=plabic_graph_from_drawing(4,colors,edges,positions)
    @test trip_permutation(G)==[-1,-2,3,4]
    @test is_reduced(G)

    original=layout_plabic_graph(plabic_graph([3,4,1,2]))
    rebuilt=plabic_graph_from_drawing(original.n,
                                      original.colors[original.n+1:end],
                                      original.edges,
                                      original.positions[original.n+1:end])
    @test trip_permutation(rebuilt)==[3,4,1,2]
    @test is_reduced(rebuilt)
    @test_throws ArgumentError plabic_graph_from_drawing(2,[:black],[(1,3)],[(0.0,0.0)])

    # A boundaryless component is removed explicitly; the reducer must not
    # replace the user's graph with a canonical bridge representative.
    extra_colors=[:white,:white,:black,:black,:black,:white]
    extra_positions=[(0.0,-.8),(-.8,0.0),(0.0,.8),(.8,0.0),(-.12,.08),(.12,.08)]
    extra_edges=[(1,5),(2,6),(3,7),(4,8),(9,10)]
    explicitly_reduced=plabic_graph_from_drawing(4,extra_colors,extra_edges,extra_positions)
    @test trip_permutation(explicitly_reduced)==[-1,-2,3,4]
    @test is_reduced(explicitly_reduced)
    @test all(all(isapprox.(actual,expected;atol=1e-12)) for (actual,expected) in
              zip(explicitly_reduced.positions[1:4],
                  [(0.0,-1.0),(-1.0,0.0),(0.0,1.0),(1.0,0.0)]))
    _,component_stages=Positroids._drawing_reduction_animation(
        4,extra_colors,extra_edges,extra_positions)
    @test any(occursin("isolated dipole",lowercase(caption))
              for (caption,_) in component_stages)
    key=Positroids._explicit_reduction_key(explicitly_reduced)
    @test key[2] == Tuple(sort(explicitly_reduced.edges))
    differently_sized=plabic_graph([3,4,1,2])
    visited=Set{Any}([key])
    push!(visited,Positroids._explicit_reduction_key(differently_sized))
    @test length(visited)==2

    # The interactive reducer retains each local operation for continuous
    # interpolation instead of returning only the final graph.
    original=layout_plabic_graph(plabic_graph([3,4,1,2]);iterations=20,restarts=1)
    animated_colors=copy(original.colors[5:end])
    animated_positions=copy(original.positions[5:end])
    animated_edges=copy(original.edges)
    boundary_edge=findfirst(e->e[1]<=4 || e[2]<=4,animated_edges)
    u,v=animated_edges[boundary_edge]
    boundary=u<=4 ? u : v
    internal=u<=4 ? v : u
    midpoint=((original.positions[boundary][1]+original.positions[internal][1])/2,
              (original.positions[boundary][2]+original.positions[internal][2])/2)
    push!(animated_colors,original.colors[internal]==:black ? :white : :black)
    push!(animated_positions,midpoint)
    inserted=4+length(animated_colors)
    animated_edges[boundary_edge]=(boundary,inserted)
    push!(animated_edges,(inserted,internal))
    animated,stages=Positroids._drawing_reduction_animation(
        4,animated_colors,animated_edges,animated_positions)
    @test trip_permutation(animated)==[3,4,1,2]
    @test is_reduced(animated)
    @test any(occursin("bivalent",lowercase(caption)) for (caption,_) in stages)

    # The two rising diagonal edges in the four-boundary-vertex example are
    # precisely the edges suppressed with their bivalent endpoints.
    picture_colors=[:black,:white,:white,:black]
    picture_positions=[(-1/3,0.0),(0.0,1/3),(0.0,-1/3),(1/3,0.0)]
    picture_edges=[(2,5),(5,6),(6,3),(5,7),(7,1),(7,8),(8,4)]
    picture_graph,picture_stages=Positroids._drawing_reduction_animation(
        4,picture_colors,picture_edges,picture_positions)
    @test length(picture_stages)==2
    @test all(occursin("Suppress one bivalent vertex",caption)
              for (caption,_) in picture_stages)
    @test picture_graph.edges==[(2,5),(5,6),(1,6),(3,5),(4,6)]

    # Regression reconstructed from the user's 10-boundary-vertex drawing.
    # Its explicit reductions change the trip permutation twice; freezing the
    # initial permutation incorrectly compares dimensions 23 and 19 forever.
    screenshot_colors,screenshot_edges,screenshot_positions=
        screenshot_reduction_fixture()
    screenshot_graph,screenshot_stages=Positroids._drawing_reduction_animation(
        10,screenshot_colors,screenshot_edges,screenshot_positions)
    @test is_reduced(screenshot_graph)
    @test trip_permutation(screenshot_graph)==[4,6,8,9,1,10,2,3,5,7]
    @test length(screenshot_stages)>=20
    @test count(stage->occursin("Complete one facial square move",stage[1]),
                screenshot_stages)==4
end

@testset "plabic force layout" begin
    G = plabic_graph([3, 4, 1, 2])
    H = layout_plabic_graph(G; iterations=50, restarts=2)
    @test H.edges == G.edges
    @test H.colors == G.colors
    @test H.positions[1:G.n] == G.positions[1:G.n]
    @test all(hypot(p...) <= 1.000001 for p in H.positions[(G.n+1):end])
    @test_throws ArgumentError layout_plabic_graph(G; restarts=0)
    @test_throws ArgumentError layout_plabic_graph(G; method=:unknown)
    G58=plabic_graph([4,5,6,7,8,1,2,3])
    H58=layout_plabic_graph(G58;method=:harmonic)
    radii=hypot.(first.(H58.positions[9:end]),last.(H58.positions[9:end]))
    @test sum(radii)/length(radii) < 0.5
    @test Positroids._edge_crossing_count(H58.positions,H58.edges)==0
end


@testset "cell parametrizations" begin
    for n in 1:5, k in 0:n, p in decorated_permutations(n, k)
        B = bridge_parametrization(p)
        N = boundary_measurement_parametrization(p)
        d = dimensionOfPermutation(k, n, p)
        @test length(parameter_names(B)) == d
        @test length(parameter_names(N)) == d
        weights = [1.1 + 0.137j for j in 1:d]
        A = parametrization_matrix(B, weights)
        C = parametrization_matrix(N, weights)
        @test A == C
        @test size(A) == (k, n)
        expected = fromDecoratedPermToPositroid(k, n, p)
        support = Set(collect(I) for (I, minor) in plucker_coordinates(A)
                      if abs(minor) > 1e-8)
        @test support == expected
        @test all(x >= -1e-8 for x in values(plucker_coordinates(A)))
    end

    B = bridge_parametrization([3, 4, 1, 2])
    @test parameter_names(B) == [:α1, :α2, :α3, :α4]
    @test size(symbolic_matrix(B)) == (2, 4)
    symbolic = plucker_coordinates(B)
    @test symbolic[(1, 2)] == "1"
    @test symbolic[(2, 3)] == "α3*α4"
    @test symbolic[(2, 4)] == "α1*α2 + α1*α4"
    N = boundary_measurement_parametrization([3, 4, 1, 2])
    @test plucker_coordinates(N)[(2, 3)] == "w3*w4"
    Z = [1 0 0; 0 1 0; 0 0 1; 1 -1 1]
    @test N * Z == ["1" "w2 + w4" "w2*w3";
                     "w1" "1 - w1" "w1 + w3"]
    @test parametrization_matrix(bridge_parametrization([-1,-2,3,4]))==
          [1.0 0.0 0.0 0.0;0.0 1.0 0.0 0.0]
    @test_throws DimensionMismatch symbolic_product(N, ones(Int, 3, 2))
    @test_throws ArgumentError parametrization_matrix(B, [1.0])
end


@testset "positroid poset children" begin
    function exhaustive_children(parent)
        n=length(parent)
        k=decorated_excedances(parent)
        d=dimensionOfPermutation(k,n,parent)
        d==0 && return Vector{Vector{Int}}()
        parent_bases=fromDecoratedPermToPositroid(k,n,parent)
        result=Vector{Vector{Int}}()
        for candidate in decorated_permutations(n,k)
            dimensionOfPermutation(k,n,candidate)==d-1 || continue
            child_bases=fromDecoratedPermToPositroid(k,n,candidate)
            isempty(setdiff(child_bases,parent_bases)) || continue
            push!(result,candidate)
        end
        sort!(result)
        return result
    end

    # The private source-convention generator agrees with the defining
    # exhaustive basis-containment calculation. The public function converts
    # target/trip permutations in and out, so callers never need `invperm`.
    for n in 1:4, k in 0:n, parent in decorated_permutations(n,k)
        source_children=exhaustive_children(parent)
        @test Positroids._immediate_children_source(parent)==source_children
        target_parent=Positroids._decorated_inverse(parent)
        expected_target=sort([Positroids._decorated_inverse(child)
                              for child in source_children])
        @test immediate_children(target_parent)==expected_target
    end

    asymmetric=[4,3,1,5,2]
    asymmetric_children=immediate_children(asymmetric)
    @test length(asymmetric_children)==7
    @test [4,-2,1,5,3] ∉ asymmetric_children
    @test all(Positroids._target_rank(child)==2 for child in asymmetric_children)
    @test Set(asymmetric_children)==Set([
        [-1,3,4,5,2],[4,1,3,5,2],[2,3,1,5,4],[4,2,1,5,3],
        [4,3,1,2,5],[5,3,1,4,2],[4,3,2,5,1]])
    @test length(boundary_cells(asymmetric))==50
    @test boundary_f_vector(asymmetric)==[8,18,17,7]

    p = [3, 4, 1, 2]
    children = immediate_children(p)
    @test !isempty(children)
    k = decorated_excedances(p)
    d = dimensionOfPermutation(k, length(p), p)
    @test all(dimensionOfPermutation(k, length(p), q) == d-1 for q in children)
    @test all(is_child(q, p) === true for q in children)
    @test immediate_children([-1, -2]) == Vector{Vector{Int}}()
    large_children=immediate_children([5,6,7,8,1,2,3,4])
    @test length(large_children)==8
    @test all(dimensionOfPermutation(4,8,q)==15 for q in large_children)

    boundary=boundary_cells(p)
    @test length(boundary)==32
    @test length(Set(Tuple.(boundary)))==length(boundary)
    @test all(is_child(q,p)===true for q in boundary)
    boundary_dimensions=[Positroids._target_cell_dimension(q) for q in boundary]
    @test issorted(boundary_dimensions)
    @test boundary_f_vector(p)==[6,12,10,4]
    @test boundary_f_vector(p;include_cell=true)==[6,12,10,4,1]
    @test boundary_cells([-1,-2])==Vector{Vector{Int}}()
    @test boundary_f_vector([-1,-2])==Int[]
    @test boundary_f_vector([-1,-2];include_cell=true)==[1]
    @test_throws ArgumentError boundary_cells(p;max_cells=3)
end

@testset "projection Jacobian ranks on positroid boundaries" begin
    # Projection routines use the user's target/anti-excedance convention by
    # default. Z is [I_4; -1 1 -1 1].
    p=[4,3,1,5,2]
    Z=[1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1;
      -1 1 -1 1]

    parent=projection_jacobian_report(p,Z;samples=1)
    @test parent.permutation==p
    @test parent.source_permutation==[3,5,2,1,4]
    @test parent.permutation_convention==:target
    @test parent.rank_k==2
    @test parent.dimension==4
    @test parent.jacobian_rank==4
    @test parent.certificate==:full_rank

    facets=boundary_projection_jacobian_report(p,Z;samples=1)
    @test length(facets)==7
    bad=filter(result->!result.full_rank,facets)
    @test length(bad)==1
    @test only(bad).permutation==[-1,3,4,5,2]
    @test only(bad).source_permutation==[-1,5,2,3,4]
    @test only(bad).jacobian_rank==2
    @test only(bad).rank_upper_bound==2
    @test only(bad).certificate==:rank_drop

    boundary=boundary_projection_jacobian_report(p,Z;strata=:boundary,samples=1)
    @test length(boundary)==50
    lower_faces=filter(result->result.dimension<3,boundary)
    @test all(result->result.full_rank,lower_faces)

    poset=projection_boundary_poset(p,Z;samples=1)
    @test length(poset.nodes)==51
    @test poset.parent==Tuple(p)
    @test poset.source_parent==(3,5,2,1,4)
    @test poset.permutation_convention==:target
    @test poset.f_vector==[8,18,17,7,1]
    @test poset.nodes[1].label=="Pi"
    @test all(edge->poset.nodes[edge.upper].dimension==
                    poset.nodes[edge.lower].dimension+1,poset.covers)
    @test length(poset.covers)==length(Set((edge.upper,edge.lower)
                                           for edge in poset.covers))
    @test only(filter(node->node.certificate==:rank_drop,poset.nodes)).label=="F1"
    @test_throws DimensionMismatch projection_jacobian_report(p,Z[1:4,:])
    @test_throws ArgumentError boundary_projection_jacobian_report(p,Z;strata=:unknown)
    @test_throws ArgumentError projection_jacobian_report(
        p,Z;permutation_convention=:unknown)

    source_parent=projection_jacobian_report(
        [3,5,2,1,4],Z;samples=1,permutation_convention=:source)
    @test source_parent.rank_k==2
    @test source_parent.dimension==parent.dimension
    @test source_parent.jacobian_rank==parent.jacobian_rank
end


@testset "plabic square moves" begin
    G = plabic_graph([3, 4, 1, 2])
    squares = square_cycles(G)
    @test squares == [[5, 6, 7, 8]]
    H = square_move(G, squares[1])
    @test is_bipartite(H)
    @test is_reduced(H)
    @test trip_permutation(H) == trip_permutation(G)
    @test H.colors[5:8] == reverse(G.colors[5:8])
    animated,stages=Positroids._square_move_animation_stages(G,squares[1])
    @test animated.edges==H.edges
    @test length(stages)>=3
    @test any(occursin("Switch",caption) for (caption,_) in stages)
    @test any(occursin("Reduced",caption) for (caption,_) in stages)
    @test occursin("\"caption\"",Positroids._interactive_animation_frame_json(stages[1]...))
    @test square_move(G).edges == H.edges
    @test_throws ArgumentError square_move(G, [1,2,3,4])
    L = square_move_by_label(G,[2,4];iterations=1000,restarts=8)
    @test is_reduced(L) && is_bipartite(L)
    @test L.colors == H.colors
    @test square_move(G;face_label=Set([2,4]),iterations=1000,restarts=8).edges == H.edges
    @test_throws ArgumentError square_move_by_label(G,[1,2,3])
    @test_throws ArgumentError square_move_by_label(G,[1,3];iterations=1000,restarts=8)
    @test_throws ArgumentError compare_plabic_graphs(G)

    # Explicit reduction search must allow state keys of different lengths:
    # trivalentizing this connected nonreduced graph adds vertices before its
    # square/parallel reduction removes them again.
    embedded=plabic_embedding(G;iterations=1000,restarts=8).graph
    search_graph=Positroids.PlabicGraph(
        4,vcat(embedded.colors,[:black,:white,:black]),
        vcat(embedded.edges,[(5,9),(9,10),(10,11),(5,11)]),
        vcat(embedded.positions,[(.10,-.48),(.25,-.48),(.25,-.33)]),
        copy(embedded.permutation),Tuple{Int,Int}[])
    searched,search_stages=Positroids._reduce_by_explicit_moves(
        search_graph;animate=true,max_depth=3,max_states=100)
    @test is_reduced(searched)
    @test length(search_stages)>=2
    @test any(occursin("square",lowercase(caption)) for (caption,_) in search_stages)
end


@testset "trips and face labels" begin
    p = [3,4,1,2]
    G = plabic_graph(p)
    E = plabic_embedding(G;iterations=1000,restarts=8)
    @test last.(graph_trips(E)) == abs.(p)
    labels = face_labels(E)
    @test length(labels) == 5
    @test all(length(label)==2 for label in labels)
    M = fromDecoratedPermToPositroid(2,4,p)
    @test all(label in M for label in labels)
    @test sort(labels) == sort([[2,4],[3,4],[2,3],[1,2],[1,4]])
    source_right=face_labels(E;convention=(:source,:right))
    target_left=face_labels(E;convention=(:target,:left))
    target_right=face_labels(E;convention=(:target,:right))
    @test all(length(label)==2 for label in source_right)
    @test target_left == [sort([p[i] for i in label]) for label in labels]
    @test target_right == [sort([p[i] for i in label]) for label in source_right]
    @test_throws ArgumentError face_labels(E;convention=(:banana,:left))

    p48 = [5,6,7,8,1,2,3,4]
    E48 = plabic_embedding(plabic_graph(p48);iterations=1500,restarts=16)
    labels48 = face_labels(E48)
    @test length(labels48) == 17
    @test all(length(label)==4 for label in labels48)
    @test last.(graph_trips(E48)) == p48
    users = Positroids._strand_lane_data(E48.trips)
    @test length(Positroids._strand_points(E48.trips[1],1,E48.graph.positions,users,0.02)) ==
          length(E48.trips[1])
    medial=Positroids._medial_strand_points(E48.trips[1],E48.graph.positions,
                                            E48.graph.colors,E48.graph.n)
    first_edge=E48.trips[1][1:2]
    midpoint=((E48.graph.positions[first_edge[1]][1]+E48.graph.positions[first_edge[2]][1])/2,
              (E48.graph.positions[first_edge[1]][2]+E48.graph.positions[first_edge[2]][2])/2)
    @test medial[2] == midpoint
    start_boundary=E48.graph.positions[E48.trips[1][1]]
    start_cross=start_boundary[1]*medial[1][2]-start_boundary[2]*medial[1][1]
    expected_start_sign=E48.graph.colors[E48.trips[1][2]]==:white ? 1 : -1
    @test sign(start_cross) == expected_start_sign
    end_boundary=E48.graph.positions[E48.trips[1][end]]
    end_cross=end_boundary[1]*medial[end][2]-end_boundary[2]*medial[end][1]
    expected_end_sign=E48.graph.colors[E48.trips[1][end-1]]==:white ? -1 : 1
    @test sign(end_cross) == expected_end_sign
    @test Positroids._face_barycenter([1,2,3,4],
        [(0.0,0.0),(2.0,0.0),(2.0,2.0),(0.0,2.0)]) == (1.0,1.0)
    large_center,large_font=Positroids._face_label_layout([1,2,3,4],
        [(0.0,0.0),(2.0,0.0),(2.0,2.0),(0.0,2.0)],"{1,2}";
        max_size=16,pixels_per_unit=300)
    small_center,small_font=Positroids._face_label_layout([1,2,3,4],
        [(0.0,0.0),(0.2,0.0),(0.2,0.2),(0.0,0.2)],"{1,2}";
        max_size=16,pixels_per_unit=300)
    @test all(isapprox.(large_center,(1.0,1.0)))
    @test all(isapprox.(small_center,(0.1,0.1)))
    @test small_font < large_font
    @test large_font≈16 atol=1e-3

    # Regression for the bottom-1 drawing: the face beside boundary arc 1--6
    # lies geometrically to the right of the strand beginning at 1.
    screenshot_n=6
    screenshot_colors=vcat(fill(:boundary,screenshot_n),
                           [:white,:black,:black,:white,:white,:black])
    screenshot_positions=[(cos(-2pi*(i-1)/screenshot_n-pi/2),
                           sin(-2pi*(i-1)/screenshot_n-pi/2))
                          for i=1:screenshot_n]
    append!(screenshot_positions,
            [(-.25,.28),(-.30,-.12),(.06,.12),(.02,-.28),(.38,.20),(.38,-.18)])
    screenshot_edges=[(3,7),(4,7),(7,8),(7,9),(2,8),(8,10),(9,10),
                      (9,11),(10,12),(11,12),(5,11),(6,12),(1,10)]
    screenshot_graph=Positroids._with_recomputed_trip_permutation(
        Positroids.PlabicGraph(screenshot_n,screenshot_colors,screenshot_edges,
                               screenshot_positions,collect(1:screenshot_n),
                               Tuple{Int,Int}[]))
    screenshot_embedding=plabic_embedding(screenshot_graph;method=:harmonic)
    screenshot_face=only(findall(f->1 in f && 6 in f,screenshot_embedding.faces))
    screenshot_label=face_labels(screenshot_embedding)[screenshot_face]
    @test screenshot_label==[5,6]
    @test 1 ∉ screenshot_label

    state=Positroids._interactive_state_json(plabic_graph([3,4,1,2]);
                                             iterations=1000,restarts=8)
    @test occursin("\"label\":[2,4]",state)
    @test occursin("\"movable\":true",state)
    @test occursin("\"strands\":[{\"source\":1",state)
    @test occursin("\"font_size\":",state)
    @test occursin("\"polygon\":[",state)
    @test occursin("\"permutation\":[3,4,1,2]",state)
    @test occursin("\"dual_edges\":[",state)
    @test occursin("\"dual_black_faces\":[",state)
    @test occursin("\"id\":1,\"label\":[2,4]",state)
    asymmetric_state=Positroids.JSON.parse(Positroids._interactive_state_json(
        plabic_graph([4,3,1,5,2]);iterations=700,restarts=5))
    @test asymmetric_state["k"]==2
    @test length(asymmetric_state["default_sources"])==2
    boundary_dual=Positroids._boundary_dual_positions(E)
    @test length(boundary_dual)==4
    @test all(isapprox(hypot(point...),1.0;atol=1e-12)
              for point in values(boundary_dual))
    @test !haskey(boundary_dual,1) # the central square face is internal
    @test Positroids._parse_interactive_permutation("[3, 4, 1, 2]")==[3,4,1,2]
    @test Positroids._parse_interactive_permutation("3 4 1 2")==[3,4,1,2]
    @test_throws ArgumentError Positroids._parse_interactive_permutation("[]")
    @test occursin("\"blank\":true",Positroids._interactive_blank_state_json())
    @test endswith(Positroids._interactive_history_json("{\"x\":1}",true,false),
                   "\"can_undo\":true,\"can_redo\":false}")
    morph_source=match(r"function drawMorph.*?function tweenFrames"s,
                       Positroids._INTERACTIVE_PLABIC_HTML).match
    correspondence_source=match(r"function morphCorrespondence.*?function drawMorph"s,
                                Positroids._INTERACTIVE_PLABIC_HTML).match
    @test !occursin("opacity",morph_source)
    @test occursin("targetA.id===targetB.id",morph_source)
    @test occursin("drawEdge(movedOld.get(a),movedOld.get(b))",morph_source)
    @test occursin("drawEdge(movedNew.get(a),movedNew.get(b))",morph_source)
    @test occursin("startsWith('Suppress one bivalent vertex')",correspondence_source)
    @test occursin("vertex.id===removed?anchorTarget",correspondence_source)
    @test occursin("Contract same-color edge",correspondence_source)
    @test !occursin("middle=",morph_source)
    @test !occursin("svg.style.opacity",Positroids._INTERACTIVE_PLABIC_HTML)
    @test !occursin("#graph{transition:opacity",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"assign-all-variables\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"view-rotation\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"indices-minus\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"indices-plus\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"swap-colors\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("fetch('/relabel'",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("fetch('/swap-colors'",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("inverseViewPoint(px,py)",Positroids._INTERACTIVE_PLABIC_HTML)
    relabel_source=match(r"async function cyclicallyRelabel.*?async function swapAllColors"s,
                         Positroids._INTERACTIVE_PLABIC_HTML).match
    @test occursin("render(await r.json())",relabel_source)
    @test !occursin("transitionState",relabel_source)
    @test occursin("edgeWeights.set(j+1,'t_'+(j+1))",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("assignAllVariablesButton.addEventListener('click',assignAllVariables)",
                   Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("assignAllVariablesButton.disabled=!!s.blank||",
                   Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("faceWeights.set(face.id,'s_'+index)",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("if(face.id===1)continue",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("faceMode?'Assign s₁,s₂,… to faces':'Assign t₁,t₂,… to edges'",
                   Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("if(weightMode.value==='face')for(const f of s.faces)",
                   Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("if(weightMode.value==='edge')s.edges.forEach",
                   Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"facets-toggle\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"facets-panel\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("fetch('/facets'",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("facetsToggleButton.addEventListener('click',toggleFacets)",
                   Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin(".facet-graph{display:block;width:100%;height:190px;pointer-events:none}",
                   Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"f-vector-toggle\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("id=\"f-vector-panel\"",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("fetch('/f-vector'",Positroids._INTERACTIVE_PLABIC_HTML)
    @test occursin("fVectorToggleButton.addEventListener('click',toggleFVector)",
                   Positroids._INTERACTIVE_PLABIC_HTML)
    @test !occursin("<svg",match(r"<section id=\"f-vector-panel\".*?</section>"s,
                                  Positroids._INTERACTIVE_PLABIC_HTML).match)
    facets_payload=Positroids.JSON.parse(Positroids._facets_payload(
        plabic_graph([3,4,1,2]);iterations=700,restarts=5))
    @test Int.(facets_payload["parent"])==[3,4,1,2]
    @test length(facets_payload["children"])==length(immediate_children([3,4,1,2]))
    @test all(length(facet["vertices"])>=facet["n"] for facet in facets_payload["children"])
    @test all(length(facet["edges"])>0 for facet in facets_payload["children"])
    @test sort([Int.(facet["permutation"]) for facet in facets_payload["children"]])==
          immediate_children([3,4,1,2])
    asymmetric_facets=Positroids.JSON.parse(Positroids._facets_payload(
        plabic_graph([4,3,1,5,2]);iterations=600,restarts=4))
    asymmetric_permutations=[Int.(facet["permutation"])
                             for facet in asymmetric_facets["children"]]
    @test [4,-2,1,5,3] ∉ asymmetric_permutations
    @test all(Positroids._target_rank(child)==2 for child in asymmetric_permutations)
    f_vector_payload=Positroids.JSON.parse(Positroids._f_vector_payload(
        plabic_graph([3,4,1,2]);max_cells=1000))
    @test f_vector_payload["dimension"]==4
    @test Int.(f_vector_payload["counts"])==[6,12,10,4]
    @test f_vector_payload["total_boundary_cells"]==32
    measurement=Positroids._plabic_boundary_measurement(plabic_graph([3,4,1,2]);
        mode=:edge,weights=Dict(1=>"a"),iterations=1000,restarts=8)
    @test measurement.sources==[1,2]
    @test size(measurement.matrix)==(2,4)
    @test measurement.matrix[:,measurement.sources]==["1" "0";"0" "1"]
    @test length(measurement.pluckers)==6
    alternate=Positroids._plabic_boundary_measurement(plabic_graph([3,4,1,2]);
        mode=:edge,weights=Dict(1=>"a",2=>"b"),sources=[1,3],
        iterations=1000,restarts=8)
    @test alternate.matrix[:,[1,3]]==["1" "0";"0" "1"]
    @test alternate.pluckers["2,4"]=="a + b"
    face_measurement=Positroids._plabic_boundary_measurement(plabic_graph([3,4,1,2]);
        mode=:face,weights=Dict(2=>"a",3=>"b",4=>"c",5=>"d"),
        iterations=1000,restarts=8)
    @test size(face_measurement.matrix)==(2,4)
    @test occursin("matrix_latex",Positroids._measurement_payload(
        plabic_graph([3,4,1,2]),"{\"mode\":\"edge\",\"sources\":[1,3],\"weights\":{}}";
        iterations=1000,restarts=8))

    p58=[4,5,6,7,8,1,2,3]
    G58=plabic_graph(p58)
    E58,labels58,movable58=Positroids._interactive_graph_state(G58;
                                                               iterations=1500,restarts=16)
    @test length(labels58)==16
    @test all(length(label)==Positroids._face_label_cardinality(p58,:left)
              for label in labels58)
    @test count(movable58)==5
    moved58=square_move_by_label(G58,labels58[findfirst(movable58)];
                                 iterations=1500,restarts=16)
    @test is_reduced(moved58) && is_bipartite(moved58)

    lollipop_graph=plabic_graph([-1,-2,3,4])
    lollipop_embedding=plabic_embedding(lollipop_graph)
    @test face_labels(lollipop_embedding)==[[1,2]]
    @test face_labels(lollipop_embedding;convention=(:source,:right))==[[3,4]]
    @test face_labels(lollipop_embedding;convention=(:target,:left))==[[1,2]]
    loop_areas=Float64[]
    for trip in lollipop_embedding.trips
        points=Positroids._medial_strand_points(trip,lollipop_embedding.graph.positions,
                                                lollipop_embedding.graph.colors,4)
        loop=points[2:end-1]
        area=sum(loop[j][1]*loop[mod1(j+1,length(loop))][2]-
                 loop[mod1(j+1,length(loop))][1]*loop[j][2]
                 for j in eachindex(loop))/2
        push!(loop_areas,area)
    end
    @test all(loop_areas[1:2] .< 0) # white loops are clockwise
    @test all(loop_areas[3:4] .> 0) # black loops are counterclockwise
end
