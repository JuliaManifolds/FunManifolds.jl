using Documenter, FunManifolds

#cd("docs")
makedocs(
    modules = [FunManifolds],
    authors = "Mateusz Baran",
    format = Documenter.HTML(
        prettyurls = false,
        assets = ["assets/plane1.png", "assets/sphere1.png"],
    ),
    sitename = "FunManifolds.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "man/getting-started.md",
            "man/examples.md",
            "man/geometry-intro.md",
            "man/manifolds.md",
        ],
        "Library" => Any["lib/public.md", "lib/internals.md"],
    ],
)
