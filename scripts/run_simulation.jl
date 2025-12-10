include("../src/PrionABM.jl")
using .PrionABM: 
    model_initiation as model_initiation, 
    compute_metrics as compute_metrics, 
    plot_dashboard as plot_dashboard, 
    plot_transmission_analysis as plot_transmission_analysis,
    seic_color as seic_color,
    step! as step!
using GLMakie
using CairoMakie
using Agents: abmexploration, abmvideo

function run_interactive(model)
    println("Starting interactive simulation...")
     params = Dict(
        :bd => 1e-4:1e-4:2e-3,
        :bi => 1e-7:1e-7:2e-3,
        :prion_survival => 0.9:0.01:0.9967,
    )
    GLMakie.activate!()

    fig, abmobs = abmexploration(model;
        agent_color = seic_color,
        agent_size = 0,
        heatarray = :V,
        heatkwargs = (colormap = :inferno, colorrange = (0, 750)),
        params = params
    )
    fig
end

function make_video(model, simulation_length::Int=520, title::String="CWD Spatial ABM", video_suffix::String="", heatmap_param::Symbol= :V, overlay_agents::Bool = false, heatmap_limit::Int=50)
    println("Generating video...")
    file_name = "results/cwd_evolution_" * video_suffix * ".mp4"
    agent_size = 0
    if overlay_agents
        agent_size = 4
    end
    abmvideo(file_name, model;
        agent_color = seic_color,
        agent_size = agent_size,
        frames = simulation_length,
        framerate = round(Int, simulation_length/10), # Always makes a 10 second video
        title = title,
        heatarray = heatmap_param,
        heatkwargs = (colormap = :inferno, colorrange = (0, heatmap_limit)),
    )
    println("Done!")
end

function generate_images(model, simulation_length::Int=520, image_suffix::String="")
    println("Running simulation for $(simulation_length/52) years ($(simulation_length) weeks)...")
    for week in 1:simulation_length
        step!(model, 1)
        if week % 52 == 0
            m = compute_metrics(model)
            println("  Year $(week÷52): Pop=$(m.final_population), " * 
                    "Prev=$(round(m.peak_prevalence*100, digits=1))%, " * 
                    "Direct=$(model.total_direct), Indirect=$(model.total_indirect)")
        end
    end

    dashboard_file_name = "results/cwd_dashboard_" * image_suffix * ".png"
    transmission_file_name = "results/cwd_transmission_analysis_" * image_suffix * ".png"

    # Generate and save plots
    println("\nGenerating dashboard...")
    fig1 = plot_dashboard(model)
    save(dashboard_file_name, fig1)

    println("Generating transmission analysis...")
    fig2 = plot_transmission_analysis(model)
    save(transmission_file_name, fig2)
end

function main()
    # The main simulation runner logic from the original file
    susceptible(x) = count(i == :S for i in x)
    exposed(x) = count(i == :E for i in x)
    infected(x) = count(i == :I for i in x)
    clinical(x) = count(i == :C for i in x)

    CairoMakie.activate!()

    println("Initializing model...")

    model = model_initiation(;
        nx=100,
        ny=100,
        n_deer=5000,
        n_infected=10,
        # Prion dynamics
        prion_survival=0.9967,    # ~4 year half-life
        shedding_rate=1.0,
        carcass_load=100.0,
        # Transmission parameters
        bd=9.2e-4,               # direct transmission coefficient
        bi=5.5e-5,               # indirect transmission coefficient
        k=10000,                    # aggregation parameter (0.01, 1, or 10000)
        ε=0.0001,                    # 0=density-dependent, 1=frequency-dependent
        transmission_radius=2,
    )

    # Various ways to run simulation and gather data, uncomment to use

    ## Run interactive simulation
    #run_interactive(model)

    ## Create video from simulation
    # Heatmap evolution of population density over 100 years with no agents overlaid 
    # make_video(model, 5200, "Population density heatmap", "popdensity", :population_grid, false, 10)

    # Heatmap evolution of prion density over 100 years with no agents overlaid
    make_video(model, 5200, "Prion density heatmap", "prionload", :V, false, 750)


    # generate_images(model, 520, "short")
    # generate_images(model, 5200, "long")

end

main()
