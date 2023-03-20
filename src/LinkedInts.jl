using Printf                    #

"""
    mutable struct LinkedInts
        next :: Vector{Int}
        prev :: Vector{Int}
        back :: Vector{Int}

        free_ptr :: Int
        free_cap :: Int
        root     :: Int

        block :: Vector{Int}
        size  :: Vector{Int}
    end

Block linked list, there are `length(next) = length(prev)` indices, `free_cap`
of which are free and the other are used. `root` gives the last one used, i.e.
starting from `root` and following `prev` until it gives zero should give the
`length(prev) - free_cap` used indices. Similarly, `free_ptr` gives the last
free index, i.e. starting from `free_ptr` and following `prev` until it gives
zero should give the `free_cap` free indices.

We always have `length(next) == length(prev) == length(back)`.
* `next[i]` is the minimum `j > i` such that `j` is used, or 0 if `i` is the
  last one used, i.e. `i = root`.
* `prev[i]` is the maximum `j < i` such that `j` is used, or 0 if `i` is the
  first one used.
* `back[i]` is such that `block[back[i]]` is the first index of the block
  in which `i` is.

* `free_cap` is the number of free slots.
* `free_ptr` last index of the `free_cap` free slots.
* `root` last used index.

We always have `length(block) == length(size)`.
* `block`: mapping from block index to first index of the block
* `size`: mapping from block index to length of the block
"""
mutable struct LinkedInts
    next :: Vector{Int}
    prev :: Vector{Int}
    back :: Vector{Int}

    free_ptr :: Int
    free_cap :: Int
    root     :: Int

    block :: Vector{Int}
    size  :: Vector{Int}
end

function LinkedInts(capacity=128)
    return LinkedInts(Int[], Int[], Int[],
                      0, 0, 0,
                      Int[], Int[])
end

allocatedlist(s::LinkedInts) = findall(s.block .> 0)
allocated(s::LinkedInts, id :: Int) = id > 0 && id <= length(s.block) && s.block[id] > 0
blocksize(s::LinkedInts, id :: Int) = s.size[id]
Base.length(s::LinkedInts) = length(s.next)
numblocks(s::LinkedInts) = length(s.block)

function Base.show(f::IO, s::LinkedInts)
    print(f,"LinkedInts(\n")
    @printf(f,"  Number of blocks: %d\n", length(s.block))
    @printf(f,"  Number of elements: %d\n", length(s.next))
    print(f,"  Blocks:\n")
    for i in 1:length(s.block)
        if s.block[i] > 0
            idxs = getindexes(s,i)
            print(f,"    #$i: $idxs\n")
        end
    end
    p = s.free_ptr
    freelst = Int[]
    while p > 0
        push!(freelst,p)
        p = s.prev[p]
    end
    print(f,"  Free: $freelst\n")

    println(f,"  free_ptr = $(s.free_ptr)")
    println(f,"  root     = $(s.root)")
    println(f,"  next     = $(s.next)")
    println(f,"  prev     = $(s.prev)")

    print(f,")")
end

"""
    ensurefree(s::LinkedInts, N :: Int)

Ensure that there are at least `N` elements free, and allocate as necessary.
"""
function ensurefree(s::LinkedInts, N :: Int)
    if s.free_cap < N
        num = N - s.free_cap

        cap = length(s.next)
        first = cap+1
        last  = cap+num

        append!(s.next, Int[i+1 for i in first:last])
        append!(s.prev, Int[i-1 for i in first:last])
        append!(s.back, zeros(Int, num))

        s.next[last] = 0
        s.prev[first] = s.free_ptr
        if s.prev[first] > 0
            s.next[s.prev[first]] = first
        end
        s.free_ptr = last
        s.free_cap += num

        return num
    else
        return 0
    end
end

function allocate_block(s::LinkedInts, N::Int, id::Integer)
    @assert(N > 0)
    ensurefree(s, N)
    # remove from free list
    ptre = s.free_ptr
    # ptre is the last index
    ptrb = ptre
    for i = 1:N-1
        s.back[ptrb] = id
        ptrb = s.prev[ptrb]
    end
    s.back[ptrb] = id
    s.block[id] = ptrb
    # ptrb is the first index

    prev = s.prev[ptrb]

    if prev > 0
        s.next[prev] = 0
    end

    s.free_ptr = s.prev[ptrb]
    s.free_cap -= N

    # insert into list `idx`
    s.prev[ptrb] = s.root
    if s.root > 0
        s.next[s.root] = ptrb
    end
    s.root = ptre

    #if ! checkconsistency(s)
    #    println("List = ",s)
    #    assert(false)
    #end

    return id
end

function create_block(s::LinkedInts, N::Int)
    @assert(N > 0)
    push!(s.size, N)
    push!(s.block, 0)
    @assert length(s.size) == length(s.block)
    return length(s.block)
end

"""
    newblock(s::LinkedInts, N :: Int)

Add a new block in list `idx`
"""
function newblock(s::LinkedInts, N::Int)::Int
    id = create_block(s, N)
    allocate_block(s, N, id)
    return id
end

"""
Move a block to the free list.
"""
function deleteblock(s::LinkedInts, id :: Int)
    if s.size[id] > 0
        ptrb = s.block[id]
        N = s.size[id]
        ptre = ptrb
        s.back[ptre] = 0
        for i in 2:N
            ptre = s.next[ptre]
            s.back[ptre] = 0
        end
        prev = s.prev[ptrb]
        next = s.next[ptre]

        # remove from list and clear the block id
        if s.root == ptre s.root = prev end
        if prev > 0 s.next[prev] = next end
        if next > 0 s.prev[next] = prev end

        s.size[id]  = 0
        s.block[id] = 0

        # add to free list
        if s.free_ptr > 0
            s.next[s.free_ptr] = ptrb
        end
        s.prev[ptrb] = s.free_ptr
        s.free_ptr = ptre
        s.next[ptre] = 0

        s.free_cap += N
    end
end

# TODO merge getoneindex and getindex
function getoneindex(s::LinkedInts, id :: Int)
    N = s.size[id]
    if N < 1
        error("No values at id")
    end

    s.block[i]
end

"""
    getindex(s::LinkedInts, id::Int)

Shortcut for `getindexes(s, id)[1]` when `s.size[id]` is 1.
"""
function getindex(s::LinkedInts, id::Int)
    @assert s.size[id] == 1
    @assert s.back[s.block[id]] == id
    return s.block[id]
end

"""
    getindexes(s::LinkedInts, id :: Int)

Return the vector of indices for the block `id`.
"""
function getindexes(s::LinkedInts, id :: Int)
    @assert length(s.next) == length(s.prev) == length(s.back)
    N = s.size[id]
    r = Vector{Int}(undef, N)
    p = s.block[id]
    for i in 1:N
        @assert s.back[p] == id
        r[i] = p
        p = s.next[p]
    end
    return r
end

function getindexes(s::LinkedInts, id::Int, target::Vector{Int}, offset::Int)
    N = s.size[id]
    p = s.block[id]
    for i in 1:N
        @assert s.back[p] == id
        target[i+offset-1] = p
        p = s.next[p]
    end
    return N
end

function getindexes(s::LinkedInts, ids::Vector{Int})
    N = sum(map(id -> s.size[id], ids))
    r = Vector{Int}(undef,N)
    offset = 1
    for id in ids
        offset += getindexes(s, id, r, offset)
    end
    return r
end

"""
Get a list if the currently free elements.
"""
function getfreeindexes(s::LinkedInts)
    N = s.free_cap
    r = Array{Int}(undef,N)
    ptr = s.free_ptr
    for i in 1:N
        @assert iszero(s.back[ptr])
        r[N-i+1] = ptr
        ptr  = s.prev[ptr]
    end
    r
end



"""
Get a list if the currently used elements.
"""
function getusedindexes(s::LinkedInts)
    N = length(s.next) - s.free_cap
    r = Array{Int}(undef,N)
    ptr = s.root
    for i in 1:N
        @assert !iszero(s.back[ptr])
        r[N-i+1] = ptr
        ptr  = s.prev[ptr]
    end
    r
end



"""
Check consistency of the internal structures.
"""
function checkconsistency(s::LinkedInts) :: Bool
    if length(s.prev) != length(s.next)
        return false
    end

    N = length(s.prev)


    if ! (all(i -> s.prev[i] == 0 || s.next[s.prev[i]] == i, 1:N) &&
          all(i -> s.next[i] == 0 || s.prev[s.next[i]] == i, 1:N))
        @assert(false)
    end

    mark = fill(false,length(s.prev))

    p = s.free_ptr
    while p != 0
        @assert iszero(s.back[ptr])
        mark[p] = true
        p = s.prev[p]
    end

    p = s.root
    while p != 0
        @assert(!mark[p])
        @assert !iszero(s.back[ptr])
        mark[p] = true
        p = s.prev[p]
    end

    if !all(mark)
        println(s)
        println(mark)
        @assert(all(mark))
    end

    return true
end
