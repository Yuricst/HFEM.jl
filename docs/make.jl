"""
Make documentation with Documenter.jl
"""

using Documenter

include(joinpath(dirname(@__FILE__), "../src/HFEM.jl"))


makedocs(
    clean = false,
    build = dirname(@__FILE__),
	modules  = [HFEM],
    format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "HFEM.jl",
    # options
    pages = [
		"Home" => "index.md",
        # "Tutorials" => Any[
        #     "Basics" => "basics.md",
        # ],
        "API" => "api.md",
		# "API" => Any[
		# 	"Core" => "api/api_core.md",
		# 	# "Problem Constructor" => "api/api_create_sft_problem.md",
		# 	# "Core Routines" => "api/api_core.md",
		# 	# "Sims-Flanagan Transcription" => "api/api_simsflanagan.md",
		# 	# "Plotting" => "api/api_plot.md",
		# ],
    ],
)