using IScatterSpectrum
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# Scattering volume for 230 MHz monostatic along the magnetic field
sv = IScatterSpectrum.ScatterVolume(230e6, 1e-6, 50000e-9, 0.0)

# Plasma, typical Ne, Te, and Ti
p  = IScatterSpectrum.Plasma(1.5e11, 2000e0, 1000e0, sv)

# Frequency array good for 230 MHz
f = 5.0:5:5000

s = [IScatterSpectrum.pwrspec(p, sv); IScatterSpectrum.pwrspec(1.0, p, sv)]
@test abs(s[2]-s[1])<1e-5
