"""
# Example

```
julia> Domain(1.0, 2.0, 3.0, 4.0)
Domain(1.0, 2.0, 3.0, 4.0)

julia> @define_glyphspec_show Domain

julia> Domain(1.0, 2.0, 3.0, 4.0)
Domain(minx=1.0, miny=2.0, maxx=3.0, maxy=4.0)
```
"""
macro define_show_with_fieldnames(T)
    quote
        function Base.show(io::IO, ::MIME"text/plain", x::TT) where {TT<:$(esc(T))}
            # Print concrete type with params, but limited to one line.
            tab = "    "
            limit = max(20, displaysize(io)[2])
            iob = IOBuffer()
            ioco = IOContext(iob, :compact=>true, io.dict...)
            s1 = repr(TT; context = ioco)
            show(ioco, MIME"text/plain"(),  s1; limit)
            s2 = String(take!(iob))[2:end - 1] # Remove quotes
            print(io, s2)
            if length(s2) > 20
                print(io, "\n$tab")
            end
            # Print field names and values, but limited to one line each
            print(io, "(")
            fields = fieldnames(TT)
            for (i, field) in enumerate(fields)
                s3 = string(field) * "=" * repr(getfield(x, field); context = ioco)
                show(ioco, MIME"text/plain"(),  s3; limit = limit  - 20)
                s4 = String(take!(iob))[2:end - 1] # Remove quotes
                print(io, s4)
                if i < length(fields)
                    print(io, ", ")
                    if length(s4) >= (limit - length(tab))
                        print(io, "\n$tab")
                    end
                end
            end
            print(io, ")")
        end
    end
end


"""
    elide_two_thirds(s::AbstractString; maxlen::Integer=100, mark::AbstractString=" … ") 
    ---> String

If `length(s) > maxlen`, return a version truncated around the 2 / 3 mark. The
result has total length `maxlen`, with the omission marked by `mark` (default: " … ").
Counts grapheme clusters for Unicode "safety".
"""
function elide_two_thirds(s::AbstractString; maxlen::Integer=100, mark::AbstractString=" … ")
    g = collect(Unicode.graphemes(s))
    n = length(g)
    length(mark) <= maxlen || throw(ArgumentError("maxlen must be ≥ length(mark)"))
    n <= maxlen && return s
    keep = maxlen - length(mark)
    left  = 2 * cld(keep, 3)   # ceil
    right = fld(keep, 3)       # floor
    string(join(@view g[1:left]), mark, right == 0 ? "" : join(@view g[n-right+1:n]))
end

# Specific to this module:
@define_show_with_fieldnames(Journey)
