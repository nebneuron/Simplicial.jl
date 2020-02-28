
# using Plotly
"""
Usage: PlotBettiCurves(Intervals::PersistenceIntervalsType,GraphDensity::Vector,maxdim)

"""
function PlotBettiCurves(Intervals::PersistenceIntervalsType,GraphDensity::Vector,maxdim::Int=3)
    Bettis = Intervals2Bettis(Intervals, length(GraphDensity), maxdim)
    Curves=(VERSION < v"0.7.0") ?  Array{Any}(maxdim) :  Array{Any}(undef,maxdim);
    
    for d=1:maxdim
        Curves[d]=Plotly.scatter(;x=GraphDensity, y=Bettis[d,:], mode="lines", name="beta_$d")
    end

    layout = [
      "title" => "Betti curves",
      "xaxis" => [
        "title" => "graph density \$\rho\$", "titlefont" => ["family" => "Courier New, monospace","size" => 18,"color" => "#7f7f7f"]
      ],
      "yaxis" => [
        "title" => "y Axis",
        "titlefont" => [
          "family" => "Courier New, monospace",
          "size" => 18,
          "color" => "#7f7f7f"
        ]
      ]
    ]

    Plotly.plot([Curves[d] for d=1:maxdim])
end
