#!/usr/bin/env python
#
# calculates the spectral reponse of a building represented as a single-degree-of-freedom (SDOF) system
# due to shaking caused by an earthquake.
#
# Note: this python script requires extra packages for
#     - matplotlib (version 2.1.2) : > pip install -U matplotlib
#     - numpy      (version 1.14.0): > pip install -U numpy
#     - obspy      (version 1.1.0) : > pip install -U obspy
from __future__ import print_function
import os,sys

import matplotlib as mpl
print("matplotlib version: ",mpl.__version__)

import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = 20, 5
plt.rcParams['lines.linewidth'] = 0.5

## numpy
import numpy as np
print("numpy version: ",np.__version__)

# do not show scipy warnings
import warnings
warnings.filterwarnings('ignore')

# uses obspy
# http://docs.obspy.org
# a python framework for seismology
import obspy
print("obspy version: ",obspy.__version__)
print("")

####################################################################################
# USER parameter

# show amplitude spectra
show_spectra = True

####################################################################################


# ## Theory
#
# To model the spectral response, we use a single-degree-of-freedom (SDOF) damped oscillator to mimick a building.
# The building has a specified mass $m$, stiffness $k$ and damping coefficient $c$. The spectral response is governed by the equation
#
# $$
# \partial_t^2 u(t) + 2 \xi \omega_0 \, \partial_t u(t) + \omega_0^2 \, u(t) = F(t)
# $$
#
# for the building's natural frequency $\omega_0 = \sqrt{\frac{k}{m}}$ and damping factor $\xi = \frac{c}{2 m \omega_0}$,
# driven by the earthquake ground motion $F(t)$. The building's natural period is $T_0 = \frac{2 \pi}{\omega_0}$.
# The earthquake forcing $F(t)$ is given by the ground motion acceleration recorded at the building site.<br>
#
# To solve this equation, we discrectize $t^n = n \Delta t$ and solve the equation by a Newmark time scheme.
# The Newmark scheme is a predictor-corrector scheme using displacement $u(t)$, velocity $\dot{u}(t)$ and acceleration $\ddot{u}(t)$.
#
# Each time step starts with a predictor step for displacement, velocity and acceleration:
# $$
# {u}^{n+1} = {u}^{n} + \Delta t \, \dot{u}^{n} + \frac{\Delta t^2}{2} \, \ddot{u}^{n} \\
# \dot{u}^{n+1} = \dot{u}^{n} + \frac{\Delta t}{2} \, \ddot{u}^{n} \\
# \ddot{u}^{n+1} = 0
# $$
#
# followed by solving the right-hand side of the SDOF equation for acceleration $\ddot{u}$:
# $$
# \ddot{u}^{n+1} = - \omega_0^2 \, u^{n+1} - 2 \xi \omega_0 \, \dot{u}^{n+1} + F(t^{n+1})
# $$
#
# and finished by a corrector step for velocity:
# $$
# \dot{u}^{n+1} = \dot{u}^{n+1} + \frac{\Delta t}{2} \, \ddot{u}^{n+1}
# $$
#
# Since the equation holds a velocity damping term, the acceleration $\ddot{u}^{n+1}$ will have to be multiplied
# by a factor $\frac{1}{M + \frac{\Delta t}{2} C}$ where mass matrix $M = 1$ and damping matrix $C = 2 \xi \omega_0$ for our case.


def download_data():
    # Myanmar 2025-03-28 Mw7.7 earthquake
    print("getting event information...")
    from obspy.clients.fdsn import Client

    event_time = obspy.UTCDateTime("2025-03-28T06:20:54")

    # data client
    # data centers: FDSN, IRIS, ORFEUS, IPGP, ETH, GFZ, RESIF, GEONET, USGS, ..
    c = Client("IRIS")

    # plus/minus 10 min
    starttime = event_time - 10 * 60
    endtime = event_time + 10 * 60
    Mw_min = 7.0

    catalog = c.get_events(starttime=starttime, endtime=endtime, minmagnitude=Mw_min)
    # or
    #catalog = c.get_events(eventid=3320897)
    print("Event catalog:")
    print(catalog)

    # see wilber
    # https://ds.iris.edu/wilber3/find_stations/11952284
    #
    # closest station: IU.CHTO
    print("getting station information...")

    # gets station informations
    t1 = event_time
    # adds 10 minutes
    t2 = t1 + 10 * 60

    # see wilber info:
    # CHTO loc 00: Streckeisen at depth 0.0
    #      loc 10: Streckeisen at depth 100.0
    #
    # traces of 00 seem to be clipped, better for loc 10
    # however, after removing and converting to acc/disp both look ok...
    # -> taking the one on the surface "00"

    net = "IU"    # network
    sta = "CHTO"  # station name
    loc = "00" #"*"
    cha = "BHZ"   # vertical component

    # see if station has some valuable records
    inventory = c.get_stations(network=net, station=sta,starttime=t1,endtime=t2)
    print(inventory)

    # dowloads seismogram as a stream
    st = c.get_waveforms(network=net, station=sta, location=loc, channel=cha,
                         starttime=t1,endtime=t2, attach_response=True)
    print(st)
    print("")

    return st

def get_processed_data_accel(st):
    ## data processing
    # trace statistics
    tr = st[0]

    #print(tr.stats)
    npts = tr.stats.npts      # number of samples
    dt = tr.stats.delta       # time step
    duration = npts * dt

    freq_Ny = 1.0/(2.0 * dt)  # Nyquist frequency
    freq_min = 1./duration    # minimal possible frequency for length of trace

    print("trace:")
    print("  dt       = ",dt)
    print("  duration = ",duration,"s"," = ",duration/60.,"min"," = ",duration/60./60.,"h")
    print("  Nyquist frequency = ",freq_Ny)
    print("  minimum frequency = ",freq_min)
    print("")

    # plotting traces
    st.plot();

    ## data processing
    # detrending & tapering
    st.detrend("linear")
    st.detrend("demean")
    st.taper(max_percentage=0.05)

    # instrument response removal
    # (unstable procedure)
    #
    # STS-2 instrument has a broadband range between 120s - 50 Hz
    # Streckeisen STS-2.5: range between 120s - 50 Hz (https://streckeisen.swiss/en/products/sts-2.5/)
    # frequency range
    f1 = 1.0/150.0  # low-frequency taper range (Hz)
    f2 = 1.0/120.0
    f3 = 50.0        # high-frequency taper range (Hz)
    f4 = 60.0

    # acceleration
    st_accel = st.copy()
    st_accel.remove_response(output="ACC",pre_filt=(f1,f2,f3,f4))

    print("data converted to acceleration:")
    print("  vertical acceleration min/max = ",st_accel[0].data.min(),"/",st_accel[0].data.max(),"(m)")
    print("")

    # plotting
    st_accel.plot()

    return st_accel

#
#
# ## Model spectral response
#
# Let us model the spectral response for this ground motion record.
#
# spectral response
def spectral_response(F,T_period,damping,deltat,show_info=False):

    # sets natural frequency
    omega_0 = 2.0 * np.pi / T_period

    # building damping coefficient
    xi = damping

    print("spectral response - natural period              = {:1.2f} (s)".format(T_period))
    print("                    building damping            = {:.1f} (%)".format(xi * 100.0))
    #print("                    natural (angular) frequency = {:1.2} (Hz)".format(omega_0))

    # equation: a + 2 xi omega v + omega^2 d = F
    #
    # we model: M a + C v + K d = F
    #           and M = 1
    #               C = 2 xi omega
    #               K = omega^2
    #
    # the building damping given as input here is directly xi.
    #
    # if we would know the building mass (m), stiffness (k) and viscous damping (c), we could find
    # the building damping factor xi = c / (2.0 * m * omega_0), and the natural frequency omega_0 = sqrt(k/m)
    #
    # as example: let estimate building mass & effective stiffness of the building with the given natural period
    #
    #   for buildings, we define stiffness k in [N/m] = [Pascal * m], and mass m in [kg]
    #   -> frequency omega = sqrt(k/m) has units sqrt([N/m]/[kg]) = sqrt([kg m s^(-2) / m] / [kg] )
    #                                                             = sqrt([1/s^2]) = [1/s] = [Hz]
    #
    #   stiffness k would related to the elastic modulus kappa by a certain area A (e.g., base area),
    #   and length L of the structure (e.g., building height):
    #       k = kappa * A/L
    #
    #   for wave speeds, we define c = sqrt( kappa /rho )
    #   with units c [m/s] and elastic modulus kappa in [N/m^2]=[Pascal] and density rho in [kg/m^3]
    #
    #
    #   steel buildings: we just make a rough estimate about the weight, stiffness a number of stories
    #                    for such a building
    #
    #                    weight   :  600 (floor) + 200 (top) tons per story
    #                    stiffness:  60,000 - 80,000 kN/m  smaller buildings are stiffer,
    #                                                      taller buildings more flexible
    #
    if show_info:
        # assumed effective building stiffness
        k = 80.0 * 10**6  # 80,000 kN/m = 80.0 * 10**6 N/m

        # estimated building mass:
        #   given the building natural frequency, we have omega_0 = sqrt( k / m )
        #                                                 with k = [N/m], m = [kg], omega = [Hz] -> mass m
        m = k / (omega_0**2)

        # estimated number of stories
        mass_per_story = 800. * 1000.0   # 800 tons
        number_of_stories = int(m / mass_per_story)

        # estimated damping coefficient c:
        #   xi = c / (2.0 * m * omega_0)     # building damping xi is unit-less, c has units [kg s^{-1}]
        c = xi * 2.0 * m * omega_0

        # info
        print("")
        print("                    example: assumed building stiffness  = {:.2f} (kN/m)".format(k/1000.0))
        print("                             -> building c               = {:.2e} (kg/s)".format(c))
        print("                                building mass            = {:.2f} (tons)".format(m/1000.0))
        print("                                estimated a {}-story building".format(number_of_stories))
        print("")

    # time loop
    DT = deltat
    NSTEP = len(F)

    u = np.zeros(NSTEP)     # displacement
    v = np.zeros(NSTEP)     # velocity
    a = np.zeros(NSTEP)     # acceleration

    # Newmark scheme for equations with damping:
    #  system: M a + C v + K d = F
    #  -> requires a mass matrix correction (M + DT/2 C) to solve for a in a Newmark scheme:
    #       a = 1/(M + DT/2 C) [ - K d - C v  + F]
    #
    #  in our case: a + 2 xi omega * v + omega^2 * d = F
    #  and therefore M = 1, C = 2 xi omega_0
    #  -> the mass "correction" factor becomes: 1/(1 + DT/2 C)
    #
    #  note that the forcing F is usually taken in the opposite direction as the ground motion acceleration
    #    F = - a_ground
    #
    # stability condition:
    #   the scheme we use is an explicit scheme and imposes a stability condition on the time step size:
    #   (see Hughes, 1987, book on Finite-element method, chapter 9, section 9.1.2)
    #     dt/T < C_critical / (2 pi)
    #   or
    #     dt   < T/(2pi) * C_critical = 1/omega * C_critical
    #  with constant C_critical~2 and period T
    if DT >= 2.0 / omega_0:
        print("  stability condition not satisfied: DT =",DT," should be < ",2.0/omega_0)
        return u,v,a

    fac_damp_scheme = 1.0 / (1.0 + DT/2.0 * 2.0 * xi * omega_0)

    for i in range(len(F)-1):
        # predictor
        u[i+1] = u[i] + DT * v[i] + DT**2/2.0 * a[i]
        v[i+1] = v[i] + DT/2.0 * a[i]
        a[i+1] = 0.0

        # solve
        a[i+1] = - omega_0**2 * u[i+1] - 2.0 * xi * omega_0 * v[i+1] - F[i+1]
        a[i+1] = fac_damp_scheme * a[i+1]

        # corrector
        v[i+1] = v[i+1] + DT/2.0 * a[i+1]

        # checks stability of the scheme
        if u[i+1] > 1.e50:
            print("modeling became unstable")
            u[:] = 0.0
            v[:] = 0.0
            a[:] = 0.0
            return u,v,a

    return u,v,a   # displacement,velocity,acceleration



def get_spectral_response(traceName=None,T_period=None,damping=0.05,SA_max_value=None,show=False):
    ## computes spectral response for trace
    global show_spectra

    print("*****************")
    print("Spectral response")
    print("*****************")
    print("")

    # get trace values
    if traceName == None:
        # download data as data stream
        st = download_data()

        # get acceleration
        st_accel = get_processed_data_accel(st)

        # selects acceleration trace
        tr = st_accel[0]

    else:
        # read from file
        if 1 == 0:
            # doesn't work...
            st = obspy.read(traceName,format='TSPAIR')
            print(st)
            print("")
            tr = st[0]
        else:
            # read data w/ numpy
            data = np.loadtxt(traceName)
            # time/trace data
            data_time = data[:,0]
            data_trace = data[:,1]
            # determine t0,deltat
            t0 = data_time[0]
            deltat = data_time[1] - data_time[0]

            # obspy trace object
            tr = obspy.Trace()

            # specfem file: **network**.**station**.**comp**.sem.ascii
            basename = os.path.basename(traceName)
            names = str.split(basename,".")
            #print("names: ",names)

            # sets trace stats
            tr.stats.network = names[0]
            tr.stats.station = names[1]
            tr.stats.channel = names[2]
            tr.stats.location = "00"
            tr.stats._format = "specfem ascii"
            tr.stats.ascii = {u'unit': u'M'}
            tr.stats.mseed = {u'dataquality': u'B'}
            # arbitrary starttime
            tr.stats.starttime = obspy.UTCDateTime(2025, 1, 1, 1, 0, 0)

            # sets trace data
            tr.data = data_trace
            #tr.stats.npts = n # automatically set
            # sets time step
            tr.stats.delta = deltat
            #tr.sampling_rate = 1.0/delta #automatically set
            # sets start time
            tr.stats.starttime += t0
            #print("start time: ",tr.stats.starttime)


    print("trace: ",tr)

    # show plots
    if show: tr.plot()

    # checks trace
    tr.verify()

    # plotting spectra
    if show_spectra:
        print("plotting amplitude spectra...")
        # vertical component
        #tr = st.select(component="Z").copy() # in case stream has more components (E,N,Z)
        tr_amp = tr.copy()
        # Nyquist frequency
        freq_Ny = 1.0/(2.0 * tr_amp.stats.delta)
        # Frequency domain
        FT = np.fft.rfft(tr_amp.data)
        # amplitude spectrum
        FT = np.abs(FT)
        # Frequency axis for plotting
        freqs = np.fft.rfftfreq(tr_amp.stats.npts, d=tr_amp.stats.delta)
        # Figure
        plt.title("amplitude spectra - raw data",size=10)
        plt.loglog(freqs, FT,color='black')
        plt.xlim(1e-2, freq_Ny)
        plt.xlabel("Frequency [Hz]",size=8)
        #plt.legend()
        # saves as JPEG file
        if traceName == None:
            name = f"{tr.stats.network}.{tr.stats.station}.{tr.stats.channel}.amplitude_spectra"
        else:
            name = f"{traceName}.amplitude_spectra"
        filename = name + ".jpg"
        plt.savefig(filename)
        print("")
        print("  amplitude spectra plotted as: ",filename)
        print("")
        if show: plt.show()

    # ground motion acceleration data
    F = tr.data.copy()

    # time step
    deltat = tr.stats.delta

    # plotting
    # time axis for plotting (in min)
    npts = tr.stats.npts      # number of samples
    dt = tr.stats.delta       # time step
    duration = npts * dt
    t_d = np.linspace(0,duration,npts) / (60.0)

    # clear figures
    plt.clf()

    ## building
    print("input building specs:")
    print("  damping        = ",damping)
    if T_period != None:
        print("  natural period = ",T_period)
    print("")

    # show specified single building response
    if T_period != None:
        # spectral response
        u,v,a = spectral_response(F,T_period,damping,deltat,show_info=True)

        print("")
        print("building response:")
        print("   displacement min /max = {} / {} (m)".format(u.min(),u.max()))
        print("   velocity     min /max = {} / {} (m/s)".format(v.min(),v.max()))
        print("   acceleration min /max = {} / {} (m/s^2)".format(a.min(),a.max()))
        print("")

        if show:
            # ground motion acceleration
            plt.plot(t_d, F, color='black', linewidth=1.0,label="ground motion acceleration")
            plt.legend()
            plt.xlabel("Time (min)")
            plt.ylabel("Ground motion acceleration")
            plt.show()

            # building displacement
            plt.plot(t_d, u, color='blue', linewidth=1.0,label="building displacement")
            plt.legend()
            plt.xlabel("Time (min)")
            plt.ylabel("building displacement (m)")
            plt.show()

            # building velocity
            plt.plot(t_d, v, color='green', linewidth=1.0,label="building velocity")
            plt.legend()
            plt.xlabel("Time (min)")
            plt.ylabel("building velocity (m/s)")
            plt.show()

            # building acceleration
            plt.plot(t_d, a, color='red', linewidth=1.0,label="building acceleration")
            plt.legend()
            plt.xlabel("Time (min)")
            plt.ylabel("building acceleration (m/s^2)")
            plt.show()

        # all done
        return

    # computes spectral response over a given period range
    # range of natural periods
    Np = 40
    DeltaTp = 0.2
    # range starting at 0.1 s
    periods = np.linspace(0.1,Np*DeltaTp,Np)

    # maximum spectral_response values at each period
    spectral_max_u = np.zeros(Np)
    spectral_max_v = np.zeros(Np)
    spectral_max_a = np.zeros(Np)
    spectral_max_a_total = np.zeros(Np)

    for i,T_period in enumerate(periods):
        # spectral response
        u,v,a = spectral_response(F,T_period,damping,deltat)

        # determines absolute maximum
        spectral_max_u[i] = np.abs(u).max()
        spectral_max_v[i] = np.abs(v).max()
        spectral_max_a[i] = np.abs(a).max()

        print("  maximum: displ = {} / veloc = {} / accel = {}".format(spectral_max_u[i],spectral_max_v[i],spectral_max_a[i]))
        print("")

        ## spectral response acceleration
        ## usually given as total acceleration measured in g
        # take total acceleration
        a_total = a + F

        # convert acceleration to unit of g (9.81 m/s^2 average gravitational acceleration)
        a_total /= 9.81

        spectral_max_a_total[i] = np.abs(a_total).max()


    # plot
    # spectral response displacement
    if show:
        plt.plot(periods, spectral_max_u, color='blue', linewidth=1.0,label="Spectral response displacement")
        plt.legend()
        plt.xlabel("Period (s)")
        plt.ylabel("Spectral response displacement (m)")
        plt.show()

    # spectral response velocity
    if show:
        plt.plot(periods, spectral_max_v, color='green', linewidth=1.0,label="Spectral response velocity")
        plt.legend()
        plt.xlabel("Period (s)")
        plt.ylabel("Spectral response velocity (m/s)")
        plt.show()

    # spectral response displacement
    if show:
        plt.plot(periods, spectral_max_a, color='red', linewidth=1.0,label="Spectral response acceleration")
        plt.legend()
        plt.xlabel("Period (s)")
        plt.ylabel("Spectral response acceleration (m/s^2)")
        plt.show()

    # total spectral response acceleration
    plt.plot(periods, spectral_max_a_total, color='black', linewidth=1.0,label="total Spectral response acceleration (g)")
    plt.legend()
    plt.xlabel("Period (s)")
    plt.ylabel("Spectral response acceleration (g)")
    plt.xticks(np.arange(0.0, max(periods)+DeltaTp, 0.5))
    # set maximum
    if SA_max_value != None:
        plt.ylim(0.0, SA_max_value)

    # saves as JPEG file
    if traceName == None:
        name = f"{tr.stats.network}.{tr.stats.station}.{tr.stats.channel}.spectral_response"
    else:
        name = f"{traceName}.spectral_response"
    filename = name + ".jpg"
    plt.savefig(filename)
    print("")
    print("  total spectral response plotted as: ",filename)
    print("")
    if show: plt.show()

    def save_spectral_response_data(filename,p,sa):
        with open(filename,"w") as f:
            f.write("# total spectral response acceleration (g)\n")
            f.write("# format: #period #maximum_spectral_acceleration (in g)\n")
            for i in range(len(p)):
                f.write(f"{p[i]} {sa[i]}\n")
        print("  data written to: ",filename)

    save_spectral_response_data(name + ".dat",periods, spectral_max_a_total)

    print("")
    print("done")
    print("")


def usage():
    print("./run_spectral_response.py accel_trace [--natural_period=val] [--damping=val] [--max_value=val] [--show]")
    print("with")
    print("  accel_trace     -  acceleration trace [e.g., DB.BANGKOK.BXZ.sema]")
    print("  natural_period  -  (optional) show single building with natural period (e.g., --natural_period=1.0)")
    print("  damping         -  (optional) specific building damping (e.g., --damping=0.01), default is a value of 0.05 (5% damping)")
    print("  max_value       -  (optional) maximum SA value for plotting total spectral response acceleration in (g) (e.g., --max_value=0.025)")
    print("  show            -  (optional) show figures (otherwise only console and plot to output files)")
    print("")
    sys.exit(1)

if __name__ == '__main__':
    # input
    traceName = None
    SA_max_value = None
    show = False
    # natural period
    T_period = None
    # by default, assuming 5% damping for the building(s)
    damping = 0.05

    if len(sys.argv) >= 2:
        # trace
        if "download" in sys.argv[1]:    # for downloading real data for the Myanmar earthquake
            pass
        else:
            traceName = sys.argv[1]

        # optional arguments
        for arg in sys.argv:
            # natural period
            if "--natural_period" in arg:
                T_period = float(arg.split('=')[1])
            # damping
            elif "--damping" in arg:
                damping = float(arg.split('=')[1])
            # set maximum SA for figures
            elif "--max_value" in arg:
                SA_max_value = float(arg.split('=')[1])
            elif "--show" in arg:
                show = True

    else:
        usage()

    # main routine
    get_spectral_response(traceName,T_period,damping,SA_max_value,show)
