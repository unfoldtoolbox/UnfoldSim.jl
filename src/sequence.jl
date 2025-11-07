"""
    rand_re(rng::AbstractRNG, machine::Automa.Machine)

Mimicks a reverse-regex, generating strings from regex instead of matching. Based on Automata.jl

# Arguments 
- `machine::Automa.Machine`: A Automa.Machine, typically output of `Automa.Compile(RE("mystring"))`

# Returns
- `result::String` : A string following the rules in `Automa.Machine`. `{}` are not supported, but e.g. `+`, `*`

# Examples
```julia-repl
julia> using Automa
julia> machine = Automa.compile(Automa.RegExp.RE("b+l+a+"))
julia> rand_re(MersenneTwister(2),machine)
"bbbbblllaaa"
```
"""
function rand_re(rng::AbstractRNG, machine::Automa.Machine)
    out = IOBuffer()
    node = machine.start

    while true
        if node.state ∈ machine.final_states
            (rand(rng) ≤ 1 / (length(node.edges) + 1)) && break
        end

        edge, node = rand(rng, node.edges)
        label = rand(rng, collect(edge.labels))
        print(out, Char(label))
    end

    return String(take!(out))
end

evaluate_sequencestring(rng, d::SequenceDesign) = evaluate_sequencestring(rng, d.sequence)


"""
    evaluate_sequencestring(rng, str::String)
    evaluate_sequencestring(rng, dS::SequenceDesign)
    evaluate_sequencestring(rng, dR::RepeatDesign)

Generate a sequence based on the reverse regex style string in `str`, `dS.sequence` or `dR.design.sequence`.

Directly converting to Automa.Compileis not possible, as we first need to match & evaluate the curly brackets. We simply detect and expand them.

# Arguments
- `str::String`: a string mimicking a regex, e.g. "b[lL]{3,4}a" should evaluate to e.g. "bLlLLa". E.g. "b+l*a{3,4}" should in principle evaluate to "bbbbaaa" or "bllllllllllllaaaa" - but right now we disallow `+` and `*` - we should revisit why exactly though.

# Returns
- `result::String` : a simulated string

# Examples
```julia-repl
julia> evaluate_sequencestring(MersenneTwister(1),"bla{3,4}")
"blaaaa"
```

See also [`rand_re`](@ref)
"""
function evaluate_sequencestring(rng, str::String)
    #match curly brackets and replace them
    @assert isnothing(findfirst("*", str)) && isnothing(findfirst("+", str)) "'Infinite' sequences are currently not supported."
    crly = collect(eachmatch(r"(\{[\d]+,[\d]+\})", str))
    for c in reverse(crly)
        m = replace(c.match, "{" => "", "}" => "")
        rep_minimum, rep_maximum = parse.(Int, split(m, ","))
        #@info str[c.offset-1]
        if str[c.offset-1] == ']'
            #@info "brackets found"
            bracket_end_idx = c.offset - 1
            bracket_start_idx = findlast("[", str[1:bracket_end_idx])[1]
            #@info bracket_end_idx,bracket_start_idx
            repeat_string = "[" * str[bracket_start_idx+1:bracket_end_idx-1] * "]"
        else
            bracket_start_idx = c.offset - 1
            bracket_end_idx = c.offset - 1
            repeat_string = string(str[c.offset-1])
        end

        replacement_string = repeat(repeat_string, rand(rng, rep_minimum:rep_maximum))
        #@info "rep" replacement_string
        str =
            str[1:bracket_start_idx-1] *
            replacement_string *
            str[bracket_end_idx+length(c.match)+1:end]
        #        @debug str
    end
    return rand_re(rng, Automa.compile(RE(str)))
end

evaluate_sequencestring(rng, d::RepeatDesign) = evaluate_sequencestring(rng, d.design)
