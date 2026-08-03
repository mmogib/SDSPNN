# Static syntax check for every Julia script. This does not include or execute them.

function main()
    files = sort(filter(path -> endswith(path, ".jl"), readdir(@__DIR__; join=true)))
    for file in files
        Meta.parseall(read(file, String); filename=file)
        println("PARSE OK ", basename(file))
    end
end

main()
