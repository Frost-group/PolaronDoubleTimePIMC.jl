using PolaronDoubleTimePIMC
using Documenter

DocMeta.setdocmeta!(PolaronDoubleTimePIMC, :DocTestSetup, :(using PolaronDoubleTimePIMC); recursive=true)

makedocs(;
    modules=[PolaronDoubleTimePIMC],
    authors="Jarvist Moore Frost",
    repo="https://github.com/Frost-group/PolaronDoubleTimePIMC.jl/blob/{commit}{path}#{line}",
    sitename="PolaronDoubleTimePIMC.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Frost-group.github.io/PolaronDoubleTimePIMC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Frost-group/PolaronDoubleTimePIMC.jl",
    devbranch="main",
)
