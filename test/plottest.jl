const defcolors = ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                   "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"]

sv = IScatterSpectrum.ScatterVolume(230e6, 1e-6, 50000e-9, 0.0)
p  = IScatterSpectrum.Plasma(1.5e11, 2000., 1000., sv)
f  = 5.0:5:5000
freq = [0; f]
tr   = GenericTrace[]
push!(tr, scatter(x=freq,
                  y=[IScatterSpectrum.pwrspec(p, sv);
                     [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                  line_color=defcolors[1],
                  name="Te=2000 K"))
p  = IScatterSpectrum.Plasma(1.5e11, 1000., 1000., sv)
push!(tr, scatter(x=freq,
                  y=[IScatterSpectrum.pwrspec(p, sv);
                     [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                  line_color=defcolors[2],
                  name="Te=1000 K"))
p  = IScatterSpectrum.Plasma(1.5e11, 1500., 1000., sv)
push!(tr, scatter(x=freq,
                  y=[IScatterSpectrum.pwrspec(p, sv);
                     [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                  line_color=defcolors[3],
                  name="Te=1500 K"))
p  = IScatterSpectrum.Plasma(1.5e11, 3000., 1000., sv)
push!(tr, scatter(x=freq,
                  y=[IScatterSpectrum.pwrspec(p, sv);
                     [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                  line_color=defcolors[4],
                  name="Te=3000 K"))
layout = Layout(xaxis_title="Frequency in Hz",
                yaxis_title="Power",
                font_size=18,
                legend_x=0.85,
                title="Collisionless Incoherent Scatter Spectrum")
plt = plot(tr, layout)
PlotlyJS.savefig(plt, "/home/scb/code/julia/pkg/v0.6/IScatterSpectrum/test/CollisionlessPlasma.pdf")
PlotlyJS.savefig(plt, "/home/scb/code/julia/pkg/v0.6/IScatterSpectrum/test/CollisionlessPlasma.png")
tr2 = [tr[k] for k in [1, 2, 4]]
sv = IScatterSpectrum.ScatterVolume(230e6, 1e-6, 50000e-9, 70*π/180)
p  = IScatterSpectrum.Plasma(1.5e11, 2000., 1000., sv)
push!(tr2, scatter(x=freq,
                   y=[IScatterSpectrum.pwrspec(p, sv);
                      [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                   line_color=defcolors[1],
                   line_dash="dash",
                   name="α=70 deg"))
p  = IScatterSpectrum.Plasma(1.5e11, 1000., 1000., sv)
push!(tr2, scatter(x=freq,
                   y=[IScatterSpectrum.pwrspec(p, sv);
                      [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                   line_color=defcolors[2],
                   line_dash="dash",
                   showlegend=false))
p  = IScatterSpectrum.Plasma(1.5e11, 3000., 1000., sv)
push!(tr2, scatter(x=freq,
                   y=[IScatterSpectrum.pwrspec(p, sv);
                      [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                   line_color=defcolors[4],
                   line_dash="dash",
                   showlegend=false))
sv = IScatterSpectrum.ScatterVolume(230e6, 1e-6, 50000e-9, 30*π/180)
p  = IScatterSpectrum.Plasma(1.5e11, 2000., 1000., sv)
push!(tr2, scatter(x=freq,
                   y=[IScatterSpectrum.pwrspec(p, sv);
                      [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                   line_color=defcolors[1],
                   line_dash="dot",
                   name="α=30 deg"))
p  = IScatterSpectrum.Plasma(1.5e11, 1000., 1000., sv)
push!(tr2, scatter(x=freq,
                   y=[IScatterSpectrum.pwrspec(p, sv);
                      [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                   line_color=defcolors[2],
                   line_dash="dot",
                   showlegend=false))
p  = IScatterSpectrum.Plasma(1.5e11, 3000., 1000., sv)
push!(tr2, scatter(x=freq,
                   y=[IScatterSpectrum.pwrspec(p, sv);
                      [IScatterSpectrum.pwrspec(freq, p, sv) for freq in f]],
                   line_color=defcolors[4],
                   line_dash="dot",
                   showlegend=false))
plt = plot(tr2, layout)
PlotlyJS.savefig(plt, "/home/scb/code/julia/pkg/v0.6/IScatterSpectrum/test/CollisionlessPlasmaMagField.pdf")
PlotlyJS.savefig(plt, "/home/scb/code/julia/pkg/v0.6/IScatterSpectrum/test/CollisionlessPlasmaMagField.png")
