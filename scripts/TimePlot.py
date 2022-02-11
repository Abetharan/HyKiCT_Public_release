import bisect
import numpy as np 
import os
import h5py
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import linalg
from scipy import constants
BOLTZMANN_CONSTANT = constants.value("Boltzmann constant")
ELECTRON_MASS = constants.value("electron mass")
PROTON_MASS = constants.value("proton mass")
ELECTRON_CHARGE = constants.value("elementary charge")
epsilon_0 = 8.85e-12
import seaborn as sns 

# plt.rcParams.update({'font.size': 40})
# plt.rcParams.update({'xtick.major.size':10})
# plt.rcParams.update({'xtick.minor.size':6})
# plt.rcParams.update({'ytick.major.size':10})
# plt.rcParams.update({'ytick.minor.size':6})
# plt.rcParams.update({'lines.markeredgewidth':2})
# plt.rcParams.update({'lines.markersize':2})
# plt.rcParams.update({'legend.fontsize':30})
# plt.rcParams.update({'legend.frameon': False})
# plt.rcParams.update({'lines.linewidth' : 3})
# plt.rcParams.update({
    # "text.usetex": True,
    # "font.family": "serif",
    # "font.serif": ["Palatino"],
# })
cmap = sns.color_palette("colorblind", as_cmap=True)
linestyles = ['-', "-.", "--", ":"]
directionary_of_HyKiCT_names = {"Te" : "ELECTRON_TEMPERATURE", "Tr" : "RAD_TEMPERATURE", "Ti" : "ION_TEMPERATURE",
                            "Rho" : "DENSITY", "Ne" : "ELECTRON_NUMBER_DENSITY", "Ni" : "ION_NUMBER_DENSITY", "Zbar" : "ZBAR", "qe" : "ELECTRON_HEAT_FLOW_X",
                        "He" : "ELECTRON_HEAT_CONDUCTION","Multi" : "qe.txt", "J" : "RAD_FREE_FREE_EMISSION", "A" : "RAD_FREE_FREE_ABSORB","InvBrem" : "INVERSE_BREM",
                        "SumRadPower" : "n/a", "Er" : "RAD_ENERGY_DENSITY", "AllRadGroup" : "RAD_MULTI_GROUP", "CoulombLog": "COULOMB_LOG","Kappa_P" : "PLANCKOPACITIES", "Kappa_R" : "ROSSLANDOPACITIES",
                        "V" : "VELOCITY", "Pt" : "TOTAL_PRESSURE"}
cgs_units  ={"x" : "x/cm", "Te" : "Te/eV", "Tr" : "Tr/eV", "Ti" : "Ti/eV", "V" : "v/ms^-1", "Pt" : "P_t / Pa",
            "Rho" : "rho/gcc", "Ne" : "ne/nc", "Ni" : "ni/cm^-3", "Zbar" : "ZBAR", "qe" : "$q_e/10^{16}Wm^{-2}$",
            "He" : "WKg^-1","Multi" : "q_vfp/q_hydro", "time" : "ns", "J" : "WCm^-3","A":"Wcm^-3", "InvBrem" : "WKg^-1",
            "SumRadPower":"J - A/WKg^-1", "Er" : "E_r / Wcm^-3", "AllRadGroup" : "E_r/Wcm^-3", "mass" : "Mass/gcm^-2","Kappa_R" : "m", "Kappa_P" : "m","PlasmaParam" : "U"}
cgs_conversions = {"x": 1e2, "Te" : BOLTZMANN_CONSTANT/ ELECTRON_CHARGE, "Tr" :  BOLTZMANN_CONSTANT/ ELECTRON_CHARGE, "Ti" :  BOLTZMANN_CONSTANT/ ELECTRON_CHARGE,
                    "V" : 1,"Pt": 1,"Rho" : 1e-3, "Ne" : 1e-0, "Ni" : 1e-6, "Zbar" : 1, "qe" : -1e-4*1e-12, "mass": 10,
                        "He" : 1,"Multi" : 1 , "time" : 1e9, "J" : 1, "A" : 1,"InvBrem" : 1, "SumRadPower": 1, "Er" : 1e-6, "AllRadGroup" : 1e-6, "Kappa_R" : 1,"Kappa_P" : 1, "PlasmaParam" : 1}
def return_ordered_files_hdf5(hdf5_file):
    all_times = hdf5_file['TIME'].keys()
    sorted_time = sorted(all_times, 
                        key = lambda x: int(os.path.splitext(x)[0].split('_')[-1]))
    indicies = [int(f.split('_')[-1]) for f in sorted_time]
    fluid_time = np.zeros(len(indicies))
    new_indiices = []
    for i,idx in enumerate(indicies):
        fluid_time[i] = np.array(hdf5_file["TIME/TIME_" + str(idx)]) 
        # new_indiices.append("TIME/TIME_" + str(idx))    
    return fluid_time, indicies

def gatherFluidData(nt,nx, path, finfo):
    var = np.zeros((nx,len(finfo)))
    hdf5_file = h5py.File(path, 'r')
    for i, (cycle, index, _) in enumerate(finfo):
        var[:, i] = np.array(hdf5_file["".join([cycle, "/Fluid_Output/ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_", str(index)])])
    return var

def myround(x, base=10):
    return base * round(x/base)

def find_index(dt, t_max, base = 1):
    index = []
    for t in t_max:
        index_dd = t / (dt)
        index.append(myround(index_dd, base))
    return(index)

def return_ordered_files(path, hdf5 = None):
    all_time_files = os.scandir(path)
    time_only = [f for f in all_time_files if not f.name.startswith('.')]
    sorted_heat_flows = sorted(time_only, 
                        key = lambda x: int(os.path.splitext(x.name)[0].split('_')[-1]))
    indicies = [int(os.path.splitext(f.name)[0].split('_')[-1]) for f in sorted_heat_flows]
    return indicies

def returnFluidTimes(path):
    fluid_time = []
    all_time = os.listdir("".join([path, "/TIME"]))
    time_only = [f for f in all_time if not f.startswith('.')]
    sorted_time = sorted(time_only, 
                        key = lambda x: int(os.path.splitext(x)[0].split('_')[-1]))
    for time_file in sorted_time: 
        time = np.loadtxt("".join([path, "/TIME/", time_file]))
        stripped_time_file = os.path.splitext(time_file)[0].split('_')[-1]
        fluid_time.append((stripped_time_file, time))
    cycle_time = time 
    return fluid_time

def return_ordered_files(path, hdf5 = None):
    all_time_files = os.scandir(path)
    time_only = [f for f in all_time_files if not f.name.startswith('.')]
    sorted_heat_flows = sorted(time_only, 
                        key = lambda x: int(os.path.splitext(x.name)[0].split('_')[-1]))
    indicies = [int(os.path.splitext(f.name)[0].split('_')[-1]) for f in sorted_heat_flows]
    return indicies

def findIndicies(times,path_f):
    if path_f.endswith("hdf5"):
        hdf5_file = h5py.File(path_f, "r")
        fluid_times, indicies_f = return_ordered_files_hdf5(hdf5_file)
        fluid_info = indicies_f
    else:
        fluid_info = returnFluidTimes(path_f)
        fluid_times = [time_info[1] for time_info in fluid_info]
        indicies_f = return_ordered_files(os.path.join(path_f, "ELECTRON_TEMPERATURE"))
    closest_f_indices = []

    for search_time in times:
        index = bisect.bisect_left(fluid_times, search_time)
        if index >= len(fluid_info):
            Warning("Left search area just appending last value ")
        else:
            closest_f_indices.append(fluid_info[index])

    return closest_f_indices

def conservationEnergy(path):
    # indicies = findIndicies(times, path)
    indicies = return_ordered_files(path + "/TIME/")
    energy = []
    times = []
    i = 0
    for index in indicies:
        
        time = np.loadtxt(os.path.join(path, 'TIME/TIME_' + str(index) + '.txt')) * 1e9
        velocity = np.loadtxt(os.path.join(path, 'VELOCITY/VELOCITY_' + str(index) + '.txt'))
        centered_v = np.array([(velocity[i+1] + velocity[i]) / 2 for i in range(len(velocity) - 1)])
        Te = np.loadtxt(os.path.join(path, 'ELECTRON_TEMPERATURE/ELECTRON_TEMPERATURE_' + str(index) + '.txt'))
        ne = np.loadtxt(os.path.join(path, 'ELECTRON_NUMBER_DENSITY/ELECTRON_NUMBER_DENSITY_' + str(index) + '.txt'))
        ni = np.loadtxt(os.path.join(path, 'ION_NUMBER_DENSITY/ION_NUMBER_DENSITY_' + str(index) + '.txt'))
        Ti = np.loadtxt(os.path.join(path, 'ION_TEMPERATURE/ION_TEMPERATURE_' + str(index) + '.txt'))
        Er = np.loadtxt(os.path.join(path, 'RAD_ENERGY_DENSITY/RAD_ENERGY_DENSITY_' + str(index) + '.txt'))
        intE = np.loadtxt(os.path.join(path, 'ELECTRON_INTERNAL_ENERGY/ELECTRON_INTERNAL_ENERGY_' + str(index) + '.txt'))
        intI = np.loadtxt(os.path.join(path, 'ION_INTERNAL_ENERGY/ION_INTERNAL_ENERGY_' + str(index) + '.txt'))
        # Er = 0
        rho = np.loadtxt(os.path.join(path, 'DENSITY/DENSITY_' + str(index) + '.txt'))
        workdone_boundary = 0
        # invbrem = np.loadtxt(os.path.join(path, 'InvBrem/InvBrem_' + str(index) + '.txt')) * rho
        if index == "0":
            Er[:] = 0
        # workdone_boundary += 1e-13 * (ne*BOLTZMANN_CONSTANT*Te + ni*BOLTZMANN_CONSTANT*Ti) *centered_v[-1]
        # if index !="0":
            # TotalEnergy = 0.5 * rho * centered_v**2 + ne*BOLTZMANN_CONSTANT*Te/(2/3) + ni*BOLTZMANN_CONSTANT*Ti/(2/3) + Er + workdone_boundary
        # else: 

            # TotalEnergy = 0.5 * rho * centered_v**2 + ne*BOLTZMANN_CONSTANT*Te/(2/3) + ni*BOLTZMANN_CONSTANT*Ti/(2/3) + Er
        
        TotalEnergy = 0.5 * rho * centered_v**2 + intE*rho + intI*rho + Er
        # print(np.sum(Er))
        # print(np.sum(ne*BOLTZMANN_CONSTANT*Te/(2/3) ))
        # print(np.sum(ni*BOLTZMANN_CONSTANT*Ti/(2/3) ))
        # print(np.sum(0.5 * rho * centered_v**2))
        # print("\n")
        # if i == 0:
        #     energy.append(sum(TotalEnergy) / sum(TotalEnergy))
        # else:
        #     energy.append(sum(TotalEnergy)/energy[i - 1])
        
        energy.append(np.sum(TotalEnergy))
        times.append(time) 
        i+=1
    print(energy)
    relative_diff = np.array([abs(energy[i + 1] - energy[i])/energy[i] for i in range(len(energy) - 1)])*100
    plt.plot(relative_diff)
    # plt.plot(times, energy)
    plt.ylabel("% relative total energy difference")
    plt.show()
def conservationMomentum(path, times, mass):
    indicies = findIndicies(times, path)
    mom = []
    times = []
    i = 0
    mass = np.cumsum(mass)
    mass = np.insert(mass, 0, 0)
    for indexs in indicies:
        index = indexs[0]
        time = np.loadtxt(os.path.join(path, 'TIME/TIME_' + str(index) + '.txt')) * 1e9
        velocity = np.loadtxt(os.path.join(path, 'VELOCITY/VELOCITY_' + str(index) + '.txt'))
        centered_v = np.array([(velocity[i+1] + velocity[i]) / 2 for i in range(len(velocity) - 1)])
        rho = np.loadtxt(os.path.join(path, 'DENSITY/DENSITY_' + str(index) + '.txt'))
        Totalmom = rho * centered_v
        # if i == 0:
        #     energy.append(sum(TotalEnergy) / sum(TotalEnergy))
        # else:
        #     energy.append(sum(TotalEnergy)/energy[i - 1])
        mom.append(np.sum(Totalmom))
        times.append(time) 
        i+=1
    relative_difff = np.array([abs(mom[i+1] - mom[i])/mom[i] for i in range(len(mom) - 1)])*100
    plt.plot(relative_difff)
    plt.ylabel("% relative total mom difference")
    plt.show()

def HyKiCT_plotter(path, var, times, name, linestyle = '-',color = "c", give_time = None, wavelength =None, ng = 1, nx =1):
    indicies = findIndicies(times, path)
    if path.endswith("hdf5"):
        hdf5 = True
        hdf5_file = h5py.File(path, "r")
    else:
        hdf5 = False
    for indexs in indicies:
        xlabel = 'x/cm'
        if hdf5:
            time = np.array(hdf5_file["TIME/TIME_" + str(indexs)])
            if var in ['qe','qi', 'v', "Multi", "spitzer"]:
                grid_x = np.array(hdf5_file["CELL_WALL_X/CELL_WALL_X_" + str(indexs)])
            else:
                grid_x = np.array(hdf5_file["CELL_CENTRE_X/CELL_CENTRE_X_" + str(indexs)])
            if var == "PlasmaParam":
                ylabel = 'u'
                fluid_Te = np.array(hdf5_file["".join([directionary_of_HyKiCT_names["Te"], "/", directionary_of_HyKiCT_names["Te"], "_" + str(indexs)])])
                fluid_ne = np.array(hdf5_file["".join([directionary_of_HyKiCT_names["Ne"], "/", directionary_of_HyKiCT_names["Ne"], "_" + str(indexs)])])
                debye = pow((epsilon_0*BOLTZMANN_CONSTANT * fluid_Te) / (fluid_ne * ELECTRON_CHARGE *ELECTRON_CHARGE), 0.5)
                Var = 1/(fluid_ne * pow(debye, 3))
            else:
                Var = np.array(hdf5_file["".join([directionary_of_HyKiCT_names[var], "/", directionary_of_HyKiCT_names[var], "_" + str(indexs)])])
        else:
            time = np.loadtxt(path + "/TIME/TIME_" + indexs[0] + ".txt")
            if var in ['qe','qi', 'v', "Multi", "spitzer"]:
                grid_x = np.loadtxt(path + "/CELL_WALL_X/CELL_WALL_X_" + indexs[0] + ".txt")
            else:
                grid_x = np.loadtxt(path + "/CELL_CENTRE_X/CELL_CENTRE_X_" + indexs[0] + ".txt")
            if var == "PlasmaParam":
                ylabel = 'u'
                fluid_Te = np.loadtxt(path + "".join(["/", directionary_of_HyKiCT_names["Te"], "/", directionary_of_HyKiCT_names["Te"], "_" + indexs[0] + ".txt"]))
                fluid_ne = np.loadtxt(path + "".join(["/", directionary_of_HyKiCT_names["Ne"], "/", directionary_of_HyKiCT_names["Ne"], "_" + indexs[0] + ".txt"]))
                debye = pow((epsilon_0*BOLTZMANN_CONSTANT * fluid_Te) / (fluid_ne * ELECTRON_CHARGE *ELECTRON_CHARGE), 0.5)
                Var = 1/(fluid_ne * pow(debye, 3))
            else:
                Var = np.loadtxt(path + "".join(["/", directionary_of_HyKiCT_names[var], "/", directionary_of_HyKiCT_names[var], "_" + indexs[0] + ".txt"]))
        if var == "Ne":
            nc = 1.1E15/pow(wavelength,2)
            Var /= nc
        # if var == "qe":
        #     fluid_Te = np.loadtxt(path + "".join(["/", directionary_of_HyKiCT_names["Te"], "/", directionary_of_HyKiCT_names["Te"], "_" + indexs[0] + ".txt"]))
        #     fluid_ne = np.loadtxt(path + "".join(["/", directionary_of_HyKiCT_names["Ne"], "/", directionary_of_HyKiCT_names["Ne"], "_" + indexs[0] + ".txt"]))
        #     boundary_Te = np.array([(fluid_Te[i] + fluid_Te[i+1])*0.5 for i in range(len(fluid_Te) - 1)])
        #     boundary_ne = np.array([(fluid_ne[i] + fluid_ne[i+1])*0.5 for i in range(len(fluid_ne) - 1)])
        #     qfs = pow(ELECTRON_MASS, -0.5) * BOLTZMANN_CONSTANT * np.sqrt(BOLTZMANN_CONSTANT)* boundary_ne * boundary_Te*np.sqrt(boundary_Te)
            # Var[1:-1] /= qfs

        # plt.plot(grid_x * cgs_conversions["x"], Var * cgs_conversions[var],linestyle, color = color, label = "{} t = {:.1f} ns".format(name, time* cgs_conversions['time']))
        # if var =="AllRadGroup":
        #     plt.semilogy(Var,linestyle, label = "{} t = {:.2f} ns".format(name, time * cgs_conversions['time']))

        if var == "AllRadGroup" or var == "Kappa_R" or var == "Kappa_P":
            if var == "Kappa_R" or var == "Kappa_P":
                scaling = np.array(hdf5_file["".join([directionary_of_HyKiCT_names["Rho"], "/", directionary_of_HyKiCT_names["Rho"], "_" + str(indexs)])])
                for i in range(ng):
                    new_var = Var[nx*i:nx*(i + 1)] * cgs_conversions[var] * scaling
                    plt.plot(1/new_var ,linestyle, label = "{} {} t = {:.2f} ns".format(name,i + 1, time * cgs_conversions['time']))
            else:
                for i in range(ng):
                    new_var = Var[nx*i:nx*(i + 1)] * cgs_conversions[var]
                    plt.semilogy(new_var ,linestyle, label = "{} {} t = {:.2f} ns".format(name,i + 1, time * cgs_conversions['time']))
        else:
            grid = grid_x*cgs_conversions["x"] 
            plt.plot(Var * cgs_conversions[var],linestyle, color=color, label = "{} t = {:.2f} ns".format(name, time * cgs_conversions['time']))

    plt.legend()
    plt.xlabel(cgs_units["x"])
    plt.ylabel(cgs_units[var])
if __name__ == "__main__":
    fig = plt.figure(1)
    fig.set_size_inches(15.5,10.5,forward=True)
    var ="InvBrem"
    # times = np.linspace(0e-9,4.5e-9, 10)
    times = [0.8e-9]

    # HyKiCT_plotter(path_1, var, times, "Spitzer", linestyle = '-', color = cmap[2], wavelength = 0.527e-6)
    HyKiCT_plotter(path_1, var, times, "Spitzer", linestyle = '-', color = cmap[0], wavelength = 0.8e-6)
    HyKiCT_plotter(path_2, var, times, "SNB", linestyle = '-', color = cmap[1], wavelength = 0.8e-6)

    # plt.plot(k* cgs_conversions["Te"],label = "txt")
    # plt.semilogy(k,label = "txt")
    # conservationEnergy(path_2)
    # conservationMomentum(path, times, mass)

    plt.grid(linestyle ='--', linewidth =0.5)

    plt.legend()
    plt.show()
