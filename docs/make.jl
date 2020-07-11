using Documenter, FunManifolds, Manifolds, ManifoldsBase

#cd("docs")
makedocs(
    modules = [FunManifolds, Manifolds, ManifoldsBase],
    authors = "Mateusz Baran",
    format = Documenter.HTML(prettyurls = false),
    sitename = "FunManifolds.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => Any["man/examples.md", "man/manifolds.md"],
        "Library" => Any["lib/public.md", "lib/internals.md"],
    ],
)
