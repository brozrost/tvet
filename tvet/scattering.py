from math import pi, sqrt, sin, cos, tan, log, exp, acos

deg = pi/180.0
B0 = 1.32
minh = 0.20
ming = -0.35
bartheta = 10.0 * deg

def f_lambert(f_L, mu_i, mu_e, alpha):
    return f_L

def f_lommel(f_L, mu_i, mu_e, alpha):
    if mu_i + mu_e > 0.0:
        return f_L / (mu_i + mu_e)
    else:
        return 0.0

def f_hapke(f_L, mu_i, mu_e, alpha):
    r"""
    Hapke law.
    Miroslav Broz (miroslav.broz@email.cz), Jan 14th 2023

    Reference: Spjuth (2009)
    Reference: cf. Hapke (1984)
    Reference: cf. Kuzminykh (2021), Physically based real-time rendering of the Moon.

    Notation:

    f_L             .. bi-directonal scattering function, Lambert, sr^-1
    f_L/(mu_i+mu_e) .. Lommel
    A_w             .. single-scattering albedo, 1
    P               .. Henyey-Greenstein function, 1-term, dark surface
    ming            .. asymmetry factor, g = 0 is isotropic
    H               .. Chandrasekhar function
    H(mu_i)H(mu_e)  .. multiple scattering, isotropic, bright surface
    1+B             .. opposition effect, shadow hiding
    B0              .. amplitude of o. e.
    minh            .. width of o. e.
    bartheta        .. mean slope, rad
    Sr              .. surface roughness
    i               .. incident angle, cos i = mu_i
    e               .. outgoing angle, cos e = mu_e
    alpha           .. phase angle, between \hat s and \hat o
    psi             .. azimuthal angle, between projections of \hat s, \hat o
    mu_i            .. directional cosine, incoming, cos(theta), 1
    mu_e            .. directional cosine, outgoing, 1
    mu_i_           .. effective mu_i, 1
    mu_e_           .. effective mu_e, 1
    eta_i           .. ditto for psi = 0, 1
    eta_e           .. ditto for psi = 0, 1
    xi              .. <cos bartheta>
    f               .. fraction of illumunation shadow hidden in visibility s. (or vice-versa)

    M(mu_i,mu_e)    .. multiple scattering, anisotropic; Eq. (2.21) NOT USED!
    1+B_C           .. opposition effect, coherent back-scattering; Eq. (2.27) NOT USED!
    K               .. porosity; Eq. (2.48) NOT USED!

    """
    global B, P, tanbartheta
    global mu_i_, mu_e_
    init_hapke(alpha)  # Ensure all globals are initialized

    if mu_i > 0.0 and mu_e > 0.0:
        A_w = 4.0*pi*f_L
        Stmp = Sr(mu_i, mu_e, alpha)
        f = f_L/(mu_i_+mu_e_) * ((1.0+B)*P + H(mu_i_,A_w)*H(mu_e_,A_w) - 1.0) * Stmp  # Eqs. (2.16), (2.31)
        f *= mu_e_/mu_e
    else:
        f = 0.0

    return f

def init_hapke(alpha):
    global B0, minh, ming, bartheta
    global B, P, tanbartheta

    B = B0/(1.0+1.0/minh*tan(alpha/2.0))                                # Eq. (2.26)
    P = (1.0-ming**2)/(1.0 + 2.0*ming*cos(alpha) + ming**2)**(3.0/2.0)  # Eq. (2.13)
    tanbartheta = tan(bartheta)

    return

def H(mu, A_w):
    """Chandrasekhar function."""

    gamma = sqrt(1.0-A_w)                                                     # Eq. (2.19)
    r0 = (1.0-gamma)/(1.0+gamma)                                              # Eq. (2.18)
    H = (1.0 - A_w*mu*(r0 + (1.0-2.0*r0*mu)/2.0 * log((1.0+mu)/mu)))**(-1.0)  # Eq. (2.17)
    return H

def Sr(mu_i, mu_e, alpha):
    """
    Surface roughness.

    Notes: If i, e < pi/2 - bartheta, then mu_i' ~ xi mu_i, mu_e' = xi mu_e.
    If i, e > pi/2 - bartheta, then the effective surface tilts towards
    the source or observer by about bartheta; except if both i, e are large
    and psi -> pi, the effective tilt goes to 0.

    """
    global B, P, tanbartheta
    global mu_i_, mu_e_

    cosi = mu_i
    cose = mu_e
    sini = sqrt(1.0-cosi**2)
    sine = sqrt(1.0-cose**2)
    tani = sini/cosi
    tane = sine/cose

    tmp = sini*sine
    if tmp != 0.0:
        cospsi = (cos(alpha)-cosi*cose)/tmp
    else:
        cospsi = 0.0

    psi = acos(min(max(cospsi, -1.0), 1.0))
    sinpsihalfsq = (sin(psi/2.0))**2

    xi = 1.0/sqrt(1.0 + pi*tanbartheta**2)
    f = exp(-2.0*tan(psi/2.0))

    E1i = exp(-2.0/(pi*tanbartheta*tani))
    E1e = exp(-2.0/(pi*tanbartheta*tane))
    E2i = exp(-1.0/(pi*(tanbartheta*tani)**2))
    E2e = exp(-1.0/(pi*(tanbartheta*tane)**2))

    eta_i = xi*(cosi + sini*tanbartheta * E2i/(2.0-E1i))  # Eq. (12.49)
    eta_e = xi*(cose + sine*tanbartheta * E2e/(2.0-E1e))  # Eq. (12.48)

    # Note: Signs in Eqs. corrected as in Hapke (1984), Eqs. (47), (48), (50), (51).

    if sini <= sine:
        K1 = cospsi*E2e + sinpsihalfsq*E2i
        K2 = 2.0 - E1e - psi/pi*E1i
        K3 = E2e - sinpsihalfsq*E2i

        mu_i_ = xi*(cosi + sini*tanbartheta * K1/K2)  # Eq. (2.36)
        mu_e_ = xi*(cose + sine*tanbartheta * K3/K2)  # Eq. (2.37)

        Sr = mu_i/eta_i * mu_e_/eta_e * xi/(1.0 - f + f*xi*mu_i/eta_i)  # Eq. (2.38)

    else:  # sini > sine
        K1 = cospsi*E2i + sinpsihalfsq*E2e
        K2 = 2.0 - E1i - psi/pi*E1e
        K3 = E2i - sinpsihalfsq*E2e

        mu_i_ = xi*(cosi + sini*tanbartheta * K3/K2)  # Eq. (2.39)
        mu_e_ = xi*(cose + sine*tanbartheta * K1/K2)  # Eq. (2.40)

        Sr = mu_i/eta_i * mu_e_/eta_e * xi/(1.0 - f + f*xi*mu_e/eta_e)  # Eq. (2.41)

    return Sr

def main():
    global B0, minh, ming, bartheta

    # Spjuth (2009), Fig. 2.10
    A_w = 0.23
    B0 = 1.32
    minh = 0.20
    ming = -0.35
    bartheta = 10.0*deg

    mu_i = 0.2
    mu_e = 0.4

    f_L = A_w/(4.0*pi)

    alpha = 30.0*deg

    init_hapke(alpha)

    f = f_hapke(f_L, mu_i, mu_e, alpha)

    print("f = ", f)

if __name__ == "__main__":
    main()