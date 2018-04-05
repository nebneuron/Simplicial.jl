using Documenter, Simplicial
makedocs(modules=[Simplicial])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/nebneuron/Simplicial.jl.git",
    julia  = "0.6",
    osname = "linux"
)
