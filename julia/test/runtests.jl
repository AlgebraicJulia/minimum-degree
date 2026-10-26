# Orders the graph of the README, and checks small random graphs.
using MDLib
using Test

@testset "MDLib" begin
    xadj = [1, 2, 5, 10, 11, 13, 16, 19, 21]
    adjncy = [2, 1, 3, 6, 2, 4, 6, 7, 8, 3, 6, 7, 2, 3, 5, 3, 5, 8, 3, 7]
    @test md(8, xadj, adjncy) == [2, 3, 8, 1, 5, 4, 7, 6]
    @test md(8, ones(Int, 8), xadj, adjncy) == [2, 3, 8, 1, 5, 4, 7, 6]
    @test md(Int32(8), Int32.(xadj), Int32.(adjncy)) == Int32[2, 3, 8, 1, 5, 4, 7, 6]
    @test md(0, [1], Int[]) == Int[]

    # every option gives a permutation
    for n in (1, 2, 10, 50), external in (false, true), dense in (10.0, 0.0, -1.0), mass in (false, true), twin in (false, true), compress in (false, true)
        adj = [Int[] for _ in 1:n]

        for u in 1:n, v in (u + 1):n
            if hash((n, u, v)) % 4 == 0
                push!(adj[u], v); push!(adj[v], u)
            end
        end

        xadj = cumsum([1; length.(adj)]); adjncy = reduce(vcat, adj; init = Int[])
        invp = md(n, rand(1:3, n), xadj, adjncy; external, dense, mass, twin, compress)
        @test sort(invp) == 1:n

        # heavy weights: the space for a reachable set is bounded by its
        # number of nodes, not by its weight
        invp = md(n, rand(1:1000, n), xadj, adjncy; external, dense, mass, twin, compress)
        @test sort(invp) == 1:n
    end

    # a path, a cycle and a wheel on 6 vertices, every weight 100
    for edges in ([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
            [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)],
            [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 3), (3, 4), (4, 5), (5, 6), (6, 2)])
        adj = [Int[] for _ in 1:6]

        for (u, v) in edges
            push!(adj[u], v); push!(adj[v], u)
        end

        xadj = cumsum([1; length.(adj)]); adjncy = reduce(vcat, adj)

        for external in (false, true), compress in (false, true)
            @test sort(md(6, fill(100, 6), xadj, adjncy; external, compress)) == 1:6
        end
    end
end
