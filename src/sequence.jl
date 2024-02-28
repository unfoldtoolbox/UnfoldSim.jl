function rand_re(rng::AbstractRNG, machine::Automa.Machine)
    out = IOBuffer()
    node = machine.start

    while true
        if node.state ∈ machine.final_states
            (rand() ≤ 1 / (length(node.edges) + 1)) && break
        end

        edge, node = rand(rng, node.edges)
        label = rand(rng, collect(edge.labels))
        print(out, Char(label))
    end

    return String(take!(out))
end

sequencestring(rng, d::SequenceDesign) = sequencestring(rng, d.sequence)
function sequencestring(rng, str::String)
    #match curly brackets and replace them
    @assert isnothing(findfirst("*", str)) && isnothing(findfirst("+", str)) "'infinite' sequences currently not supported"
    crly = collect(eachmatch(r"(\{[\d],[\d]\})", str))
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
        #@info str
    end
    return rand_re(rng, Automa.compile(RE(str)))
end