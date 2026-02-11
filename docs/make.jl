using MzvHoffmanAlgebra
using Documenter

DocMeta.setdocmeta!(MzvHoffmanAlgebra, :DocTestSetup, :(using MzvHoffmanAlgebra); recursive=true)

makedocs(;
    modules=[MzvHoffmanAlgebra],
    authors="いーな <kankitsumikan0424@gmail.com>, Y.K.495 <yuttarikanata@gmail.com>",
    sitename="MzvHoffmanAlgebra.jl",
    format=Documenter.HTML(;
        canonical="https://YuttariKanata.github.io/MzvHoffmanAlgebra.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/YuttariKanata/MzvHoffmanAlgebra.jl",
    devbranch="master",
)
