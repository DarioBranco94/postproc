#
# Copyright (c) 2019-2020 by University of Campania "Luigi Vanvitelli".
# Developers and maintainers: Salvatore Venticinque, Dario Branco.
# This file is part of GreenCharge
# (see https://www.greencharge2020.eu/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import pandas
import scipy.interpolate as inter
from scipy.integrate import simps





######################################################################
# Create a Database Table to store all events/devices of a scenario  #
######################################################################
def createTable():
    workingdir = config.Configuration.parameters['workingdir']
    db_filename = workingdir + '/xml/input.db'
    conn = sqlite3.connect(db_filename)
    schema_filename = '../xml/schema.sql'

    with open(schema_filename, 'rt') as f:
        schema = f.read()
    conn.executescript(schema)
    f = open(workingdir + "/xml/neighborhood.xml", "r")
    fileIntero = f.read()
    root = ET.fromstring(fileIntero)

    for house in root.findall('house'):  # READ XML FILE
        houseId = house.get('id')
        query = "insert into housecs (id) values(?)"
        cursor = conn.cursor()
        cursor.execute(query, (houseId,))

        for user in house.findall('user'):
            userId = user.get('id')
            for deviceElement in user.findall('device'):
                deviceId = deviceElement.find('id').text
                type = deviceElement.find('type').text
                name = deviceElement.find('name').text
                query = "insert into device(id,id_house,name,type,class) values(?,?,?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, houseId, name, type, type,))

            for OtherElement in user.findall('heatercooler'):
                deviceId = OtherElement.find('id').text
                name = OtherElement.find('name').text
                query = "insert into device(id,id_house,name,type,class) values(?,?,?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, houseId, name, "heatercooler", "N-S Consumer",))

            for OtherElement in user.findall('backgroundload'):
                deviceId = OtherElement.find('id').text
                name = OtherElement.find('name').text
                query = "insert into device(id,id_house,name,type,class) values(?,?,?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, houseId, name, "backgroundLoad", "N-S Consumer",))

            for OtherElement in user.findall('ecar'):
                deviceId = OtherElement.find('id').text
                name = OtherElement.find('name').text

                capacity = OtherElement.find('capacity').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "Capacity", str(capacity),))

                maxchpowac = OtherElement.find('maxchpowac').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "max_ch_pow_ac", str(maxchpowac),))

                maxchpowcc = OtherElement.find('maxchpowcc').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "max_ch_pow_cc", str(maxchpowcc),))

                maxdispow = OtherElement.find('maxdispow').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "max_dis_pow", str(maxdispow),))

                maxallen = OtherElement.find('maxallen').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "max_all_en", str(maxallen),))

                minallen = OtherElement.find('minallen').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "min_all_en", str(minallen),))

                sbch = OtherElement.find('sbch').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "sbch", str(sbch),))

                sbdis = OtherElement.find('sbdis').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "sbdis", str(sbdis),))

                cheff = OtherElement.find('cheff').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "cheff", str(cheff),))

                dis_eff = OtherElement.find('dis_eff').text
                query = "insert into staticParameter(iddevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "dis_eff", str(dis_eff),))

                query = "insert into device(id,id_house,name,type,class) values(?,?,?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, houseId, name, "ecar", "Prosumer",))

    f = open(workingdir + "/xml/loads.xml", "r")
    fileIntero = f.read();
    root = ET.fromstring(fileIntero)
    for house in root.findall('house'):
        houseId = house.get('id')
        for user in house.findall('user'):
            userId = user.get('id')
            for device in user.findall('device'):
                deviceId = device.find('id').text

                est = device.find('est').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "est", str(est),))

                lst = device.find('lst').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "lst", str(lst),))
                creation_time = device.find('creation_time').text
                type2 = device.find("type").text

                for c in Entities.listDevice:
                    if deviceId == c.id and houseId == c.house:
                        if c.type == "Consumer":
                            query = "insert into event(creation_time,idDevice,type) values(?,?,?)"
                            cursor = conn.cursor()
                            cursor.execute(query, (creation_time, deviceId, "Load",))
                        elif c.type == "Producer":
                            query = "insert into event(creation_time,idDevice,type) values(?,?,?)"
                            cursor = conn.cursor()
                            cursor.execute(query, (creation_time, deviceId, "Create PV",))

            for device in user.findall('backgroundload'):
                deviceId = device.find('id').text
                query = "insert into event(creation_time,idDevice,type) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, ("0", deviceId, "Create BG",))

            for device in user.findall('ecar'):
                deviceId = device.find('id').text
                pat = device.find('pat').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "pat", str(pat),))
                pdt = device.find('pdt').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "pdt", str(pdt),))
                aat = device.find('aat').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "aat", str(aat),))
                adt = device.find('adt').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "adt", str(adt),))
                creation_time = device.find('creation_time').text
                soc = device.find('soc').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "soc", str(soc),))
                targetSoc = device.find('targetSoc').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "targetSoc", str(targetSoc),))
                V2G = device.find('V2G').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "V2G", str(V2G),))
                priority = device.find('priority').text
                query = "insert into dinamicParameter(idDevice,key,val) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (deviceId, "targpriorityetSoc", str(priority),))
                query = "insert into event(creation_time,idDevice,type) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, (creation_time, deviceId, "ECAR event",))

            for device in user.findall('heatercooler'):
                deviceId = device.find('id').text
                query = "insert into event(creation_time,idDevice,type) values(?,?,?)"
                cursor = conn.cursor()
                cursor.execute(query, ("0", deviceId, "Create HC",))
    conn.commit()
    conn.close()


















class Intersection:
    """This class compute the intersection between two curve using different
    methods.
    """
    @classmethod
    def solve(cls, f, x):
        """
        Args:
            f:
            x:
        """
        s = np.sign(f)
        z = np.where(s == 0)[0]
        if len(z) > 0:
            return z
        else:
            s = s[0:-1] + s[1:]
            z = np.where(s == 0)[0]
            return z

    @classmethod
    def interp(cls, f, x, z):
        """
        Args:
            f:
            x:
            z:
        """
        m = (f[z + 1] - f[z]) / (x[z + 1] - x[z])
        return x[z] - f[z] / m

    @classmethod
    def intersect1(cls, x, y1,y2):
        """
        Args:
            x:
            y1:
            y2:
        """
        f = y1 - y2
        z = cls.solve(f, x)
        ans = cls.interp(f, x, z)
        return ans

    @staticmethod
    def intersect2(f1, f2, xx):
        """
        Args:
            f1:
            f2:
            xx:
        """
        def diff_func(x):
            return f1(x) - f2(x)

        roots = np.argwhere(np.diff(np.sign(diff_func(xx)))).flatten()
        return roots


class EnergyOutput:
    cons_series = {'EV', 'SH', 'HC', 'BG'}
    prod_series = {'PV'}
    sample_time = None
    self_consumption = None
    tot_consumption = None
    tot_production = None
    productions = None
    consumptions = None
    interval = 600
    folder = None



    def __init__(self, folder, interval=600):
        """
        Args:
            folder:
            interval:
        """
        self.interval = interval
        self.folder = folder + "/output"
        self.min_time = 0
        self.max_time = 0
        pass

    def load(self, rel=False):
        """
        Args:
            rel:
        """
        self.productions = self.load_series(self.prod_series, rel)
        self.consumptions = self.load_series(self.cons_series, rel)
        self.sample_time = np.arange(300, 23.5*3600, self.interval)

    def compute_production(self):
        temp_series = None
        if len(self.productions.keys()) > 0:
            for group_key in self.productions.keys():
                for ts_key in self.productions[group_key]:
                    temp_series = self.__sum_cts(temp_series, self.productions[group_key][ts_key])
            self.tot_production = temp_series
        else:
            self.tot_production = np.zeros(len(self.sample_time))

    def compute_consumption(self):
        temp_series = None
        if len(self.consumptions.keys()) > 0:
            for group_key in self.consumptions.keys():
                for ts_key in self.consumptions[group_key]:
                    temp_series = self.__sum_cts(temp_series, self.consumptions[group_key][ts_key])
            self.tot_consumption = temp_series
        else:
            self.tot_consumption = np.zeros(len(self.sample_time))

    def compute_self(self):
        self.self_consumption = np.zeros(len(self.sample_time))
        if self.tot_production is not None and self.tot_production[-1] > 0:
            for i in range(1, len(self.sample_time)):
                cons = self.tot_consumption[i] - self.tot_consumption[i-1]
                prod = self.tot_production[i] - self.tot_production[i-1]
                self_incr = 0
                if cons >= 0 and prod > 0:
                    self_incr = min(prod, cons)
                self.self_consumption[i] = self.self_consumption[i - 1] + self_incr


    def res_power(self):
        pow_cons = ce2p(self.sample_time, self.tot_consumption)
        pow_prod = ce2p(self.sample_time, self.tot_production)
        res_pow = pow_cons - pow_prod
        for i in range(len(res_pow)):
            if res_pow[i] < 0:
                res_pow[i] = 0
        return res_pow

    def __sum_cts(self, ts_sum, ts_new):
        """
        Args:
            ts_sum:
            ts_new:
        """
        return self.__sum_sub_cts(ts_sum, ts_new)

    def __sub_cts(self, ts_sum, ts_new):
        """
        Args:
            ts_sum:
            ts_new:
        """
        return self.__sum_sub_cts(ts_sum, ts_new, -1)

    # mult should be -1 to subtract, neg should be True if cumulative timeseries can decrease
    def __sum_sub_cts(self, ts_sum, ts_new, mult=1, neg=True):
        """
        Args:
            ts_sum:
            ts_new:
            mult:
            neg:
        """
        if ts_sum is None:
            ts_sum = np.zeros(len(self.sample_time))
        start_time = ts_new[0, 0]
        start_charge =ts_new[0, 1]
        xx = ts_new[:, 0] - start_time
        yy = ts_new[:, 1] - start_charge
        try:
            # ts_spline = inter.InterpolatedUnivariateSpline(xx, ts_new[:, 1])
            ts_spline = inter.interp1d(xx, yy)
        except:
            ts_spline = inter.interp1d(xx, yy)
        last_incr = 0
        for i in range(len(self.sample_time)):
            if ts_new[0, 0] <= self.sample_time[i] <= ts_new[-1, 0]:
                incr = ts_spline(self.sample_time[i]-ts_new[0, 0]) * mult
                if incr > 0 or neg:
                    last_incr = incr
                    ts_sum[i] += ts_spline(self.sample_time[i]-ts_new[0, 0])
            elif self.sample_time[i] > ts_new[-1, 0]:
                ts_sum[i] += last_incr

        return ts_sum

    def load_series(self, typ, rel=False):
        """
        Args:
            typ:
            rel:
        """
        ts_groups = {}
        for ser_typ in typ:
            series = glob.glob(self.folder + "/" + ser_typ + "/*.csv")
            ts_groups[ser_typ] = {}
            for ts_file in series:
                ts = np.genfromtxt(ts_file, delimiter=' ')
                ts[:, 0] += 7250
                ts = np.vstack((ts, ts[-1]))
                ts[-1,0] += 60
                if rel:
                    if(ser_typ == 'EV'):
                        start_time = ts[0, 0] - np.mod(ts[0, 0], 86400)
                    else:
                        start_time = ts[0, 0] - np.mod(ts[0, 0], 86400)
                    ts[:, 0] = ts[:, 0] - start_time
                if ts[0, 0] < self.min_time or self.min_time == 0:
                    self.min_time = ts[0, 0]
                if ts[-1, 0] > self.max_time:
                    self.max_time = ts[-1, 0]
                ts_groups[ser_typ][os.path.basename(ts_file)[:-4]] = ts
        return ts_groups


def ce2p(xx, yy, positive=True):
    """
    Args:
        xx:
        yy:
        positive:
    """
    yy1 = np.zeros(len(yy))
    for i in range(1, len(yy)):
        if yy[i] >= yy[i - 1] or not positive:
            yy1[i-1] = 3600 * (yy[i] - yy[i - 1]) / (xx[i] - xx[i - 1])

        else:
            yy1[i] = yy1[i - 1]

    return yy1

def ce2e(yy):
    """
    Args:
        yy:
    """
    yy1 = np.zeros(len(yy))
    for i in range(0, len(yy)-1):
        yy1[i] = yy[i+1] - yy[i]
    return yy1

def p2ce(xx, yy):
    """
    Args:
        xx:
        yy:
    """
    yy1 = np.zeros(len(yy))
    for i in range(1, len(yy)):
            yy1[i] = yy1[i-1] + yy[i] * (xx[i] - xx[i - 1])/3600
    return yy1


def plot_output(sim_output, path="."):

    """
    Args:
        sim_output:
        path:
    """
    plt.figure()
    plot_power(sim_output.productions, 'b')
    plot_power(sim_output.consumptions)
    # plt.legend()
    plt.ylabel("power (W)")
    plt.xlabel("hour")
    xlim = np.arange(sim_output.min_time, 60 * 60 * 24, 60 * 60 * 3)
    plt.xticks(xlim, [str(n).zfill(2) + ':00' for n in np.arange(int(sim_output.min_time / 3600), 24, 3)])
    plt.savefig(path+"/output_power.png")

    plt.figure()
    if sim_output.tot_production[-1] > 0:
        plt.plot(sim_output.sample_time, ce2p(sim_output.sample_time, sim_output.tot_production), 'b', linestyle='-', label="tot_production")

    if sim_output.tot_consumption[-1] > 0:
        plt.plot(sim_output.sample_time, ce2p(sim_output.sample_time, sim_output.tot_consumption), 'r', linestyle='-', label="tot_consumption")

    self_power_temp = ce2p(sim_output.sample_time, sim_output.self_consumption)
    self_power_temp[-1] = 0

    if sim_output.self_consumption[-1] > 0:
        plt.fill(sim_output.sample_time, self_power_temp, 'g', label="self_consumption")
    plt.legend()
    plt.ylabel("power (W)")
    plt.xlabel("hour")
    xlim = np.arange(sim_output.sample_time[0], 60 * 60 * 24, 60 * 60 * 3)
    plt.xticks(xlim, [str(n).zfill(2) + ':00' for n in np.arange(int(sim_output.sample_time[0] / 3600), 24, 3)])
    plt.savefig(path+"/output_self.png")

    plt.figure()
    plt.plot(sim_output.sample_time, sim_output.res_power(), 'r', label="grid_in")
    plt.legend()
    plt.ylabel("power (W)")
    plt.xlabel("hour")
    xlim = np.arange(sim_output.sample_time[0], 60 * 60 * 24, 60 * 60 * 3)
    plt.xticks(xlim, [str(n).zfill(2) + ':00' for n in np.arange(int(sim_output.sample_time[0] / 3600), 24, 3)])
    plt.savefig(path+"/grid_in.png")



def plot_power(groups, colr=None):
    """
    Args:
        groups:
        colr:
    """
    for group_key in groups.keys():
        for ts_key in groups[group_key].keys():
            ts = groups[group_key][ts_key]
            xx = ts[:, 0]
            yy = ts[:, 1]
            yy1 = np.zeros(len(yy))


            yy1 = ce2p(xx, yy)
            if yy1[0] != 0:
                xx = np.insert(xx, 0, xx[0] - 60)
                yy1 = np.insert(yy1,0 , 0)
            if colr is not None:
                plt.plot(xx, yy1, colr, linestyle='-',   label=group_key + "_" + ts_key)
            else:
                plt.plot(xx, yy1, linestyle='-', label=group_key + "_" + ts_key)

class Performance:


    @staticmethod
    def self_consumption(self_consumption, production):
        """
        Args:
            self_consumption:
            production:
        """
        if production[-1] > 0:
            return self_consumption[-1]/production[-1]
        else:
            return 0

    @staticmethod
    def average(groups):
        """
        Args:
            groups:
        """
        energy = 0
        min_time = None
        max_time = None
        for group_key in groups.keys():
            for ts_key in groups[group_key].keys():
                serie = groups[group_key][ts_key]
                energy += serie[-1, 1]
                if min_time is None or min_time > serie[0, 0]:
                    min_time = serie[0, 0]
                if max_time is None or max_time < serie[-1, 0]:
                    max_time = serie[-1, 0]

        average = 3600 * energy/(max_time-min_time)
        return average

    @staticmethod
    def peak2average(consumption):
        """
        Args:
            consumption:
        """
        avg_pow = 3600 * consumption[-1, 1]/(consumption[-1, 0] - consumption[0, 0])
        pow_serie = ce2p(consumption[:, 0], consumption[:, 1])
        return np.max(pow_serie)/avg_pow


def shift_load(shift_time, infile, outfile):
    """This code shift a timeseries :param shift_time: :param infile: :param
    outfile: :return:

    Args:
        shift_time:
        infile:
        outfile:
    """
    series = np.genfromtxt(infile, delimiter=' ')
    start_time = series[0, 0]
    series[:, 0] = series[:, 0] - start_time + shift_time
    np.savetxt(outfile, series, delimiter=' ', fmt="%d %f")


def compute_area(ffunc, xx, a, b):
    """This code compute the area below the curve ffunc, between a and b

    Args:
        ffunc:
        xx:
        a:
        b:
    """
    return simps(ffunc(xx[a:b]), xx[a:b])/3600


def ev2maxself(ev_ce, xx, res_energy, filename):
    """This function change the energy profile of an ev to optimize the usage of
    residual green energy :param ev_ce: :param xx: :param res_energy: :param
    filename: :return:

    Args:
        ev_ce:
        xx:
        res_energy:
        filename:
    """
    ev_energy = ev_ce[-1, 1]
    ev_max_en = (ev_ce[-1, 1] / (ev_ce[-1, 0] - ev_ce[0, 0])) * 600
    print("ev_demand", ev_energy)
    for i in range(len(res_energy)):
        if res_energy[i] < 0:
            res_energy[i] = 0

    charged_energy = [0]
    charged_times = []
    for i in range(len(xx)):
        if xx[i] >= 10 * 3600:
            if res_energy[i] < ev_max_en:
                usable_energy = res_energy[i]
            else:
                usable_energy = ev_max_en
            ev_energy -= usable_energy

            charged_energy.append(charged_energy[-1] + usable_energy)
            charged_times.append(xx[i])
            if ev_energy <= 0:
                print('charged at:', xx[i] / 3600, ":", xx[i] % 3600)
                break
    print("ev_residual", ev_energy)

    charged_ts = np.vstack((charged_times, charged_energy[:-1])).T
    # print(charged_ts)
    np.savetxt(filename, charged_ts, delimiter=' ', fmt="%d %f")



def callExternal(folder, out_path='.'):
    """
    Args:
        folder:
        out_path:
    """
    print(folder)

    sim_output = EnergyOutput(folder,150)
    sim_output.load(True)
    sim_output.compute_production()
    sim_output.compute_consumption()
    sim_output.compute_self()
    tot_consumption = np.vstack((sim_output.sample_time, sim_output.tot_consumption)).T
    print('tot_cons PAR:', Performance.peak2average(tot_consumption))
    print('cons PEAK:', np.max(ce2p(sim_output.sample_time, sim_output.tot_consumption)))
    res_energy = np.vstack((sim_output.sample_time, p2ce(sim_output.sample_time, sim_output.res_power()))).T
    print('res_energy PAR:', Performance.peak2average(res_energy))
    print('res_energy PEAK:',np.max(sim_output.res_power()))
    print('self consumption:', Performance.self_consumption(sim_output.self_consumption, sim_output.tot_production))
    plot_output(sim_output, out_path)






if __name__ == "__main__":
    folder = "/home/salvatore/projects/gcsimulator/docker/users/demo/Simulations/trivial/Results/12_12_15_82"
    folder = "/home/salvatore/projects/gcsimulator/docker/users/demo/Simulations/demo2"

    # folder = "/home/salvatore/projects/gcsimulator/docker/users/demo/Simulations/trivialEV/Results/12_12_15_1"
    # folder = "/home/salvatore/projects/gcsimulator/docker/users/demo/Simulations/trivial/Results/12_12_15_1"
    # shift_load(10 * 3600, folder + "/input/EV/4_run_3_1_ecar.csv",  folder + "/output/EV/4_run_3_1_ecar.csv")
    # shift_load(15 * 3600, folder + "/input/SH/10_run_1_1_dw.csv", folder + "/output/SH/10_run_1_1_dw.csv")
    # shift_load(10 * 3600, folder + "/input/SH/10_run_2_1_wm.csv", folder + "/output/SH/10_run_2_1_wm.csv")
    # shift_load(9.5 * 3600, folder + "/input/EV/10_54_.csv", folder + "/output/EV/10_54_.csv")


    sim_output = EnergyOutput(folder,150)

    sim_output.load(True)
    sim_output.compute_production()
    sim_output.compute_consumption()
    sim_output.compute_self()
    plot_output(sim_output)
    tot_consumption = np.vstack((sim_output.sample_time, sim_output.tot_consumption)).T
    tot_production = np.vstack((sim_output.sample_time, sim_output.tot_production)).T
    self_cons = np.vstack((sim_output.sample_time, sim_output.self_consumption)).T

    print('tot_cons PAR:', Performance.peak2average(tot_consumption))
    print('cons PEAK:', np.max(ce2p(sim_output.sample_time, sim_output.tot_consumption)))
    res_energy = np.vstack((sim_output.sample_time, p2ce(sim_output.sample_time, sim_output.res_power()))).T
    print('res_energy PAR:', Performance.peak2average(res_energy))
    print('res_energy PEAK:',np.max(sim_output.res_power()))
    print('self consumption:', Performance.self_consumption(sim_output.self_consumption, sim_output.tot_production))
    np.savetxt("consumption.csv", tot_consumption, delimiter=" ", fmt="%i %f")
    np.savetxt("production.csv", tot_production, delimiter=" ", fmt="%i %f")
    np.savetxt("self_consumption.csv", self_cons, delimiter=" ", fmt="%i %f")
    np.savetxt("grid_in.csv", res_energy, delimiter=" ", fmt="%i %f")


    '''
    Find the when the total consumed power intersect the maximum power (6000W)
    '''
    threshold = np.zeros(len(sim_output.sample_time)) + 6000
    tot_power = ce2p(sim_output.sample_time, sim_output.tot_consumption)
    print("Intersect", Intersection.intersect1(sim_output.sample_time, threshold, tot_power))

    func1 = inter.interp1d(sim_output.sample_time,tot_power)

    def func2(x):
        """
        Args:
            x:
        """
        return threshold
    result = Intersection.intersect2(func1, func2, sim_output.sample_time)
    print("Intersect", result)
    print(sim_output.sample_time[result])
    exit(0)

    """
    This code is used to exploit charge flexibility in order to keep the total power below a threshold
    But it is not general and must be changed.
    """
    folder = "/home/salvatore/projects/gcsimulator/docker/users/demo/Simulations/demo2"
    sim_output = EnergyOutput(folder, 150)
    sim_output.cons_series = {'HC', 'SH'}
    sim_output.prod_series = {'EV'}


    sim_output.load(True)
    sim_output.compute_production()
    sim_output.compute_consumption()
    sim_output.compute_self()

    total_power = ce2p(sim_output.sample_time, sim_output.tot_consumption)
    total_available = np.zeros((len(total_power))) +6000
    total_available -= total_power
    total_evpower = ce2p(sim_output.sample_time, sim_output.tot_production)

    ev3 = np.zeros((len(total_power)))
    ev1= np.zeros((len(total_power)))
    ev2 = np.zeros((len(total_power)))
    ev3_energy=9816
    ev2_energy=8970
    ev1_energy=8231

    for i in range(1, len(sim_output.sample_time)):
        av_pow = total_available[i] - total_power[i]
        if sim_output.sample_time[i] > 25200 and ev3_energy>0 and av_pow>0:
            ev3[i]=min(av_pow, 1621)
            av_pow-=ev3[i]
            ev3_energy-= 0.5*(ev3[i]+ev3[i-1])*(sim_output.sample_time[i]-sim_output.sample_time[i-1])/3600
            if ev3_energy<0:
                print('ev3 finish at', sim_output.sample_time[i] )
        if sim_output.sample_time[i] > 34200 and ev2_energy>0 and av_pow>0:
            ev2[i]=min(av_pow, 1546)
            av_pow-=ev2[i]
            ev2_energy-=0.5 * (ev2[i] + ev2[i - 1]) * (sim_output.sample_time[i] - sim_output.sample_time[i - 1]) / 3600
            if ev2_energy < 0:
                print('ev2 finish at', sim_output.sample_time[i])

        if sim_output.sample_time[i] > 36200 and ev1_energy>0 and av_pow>0:
            ev1[i]=min(av_pow, 3228)
            av_pow-=ev1[i]
            ev1_energy-=0.5 * (ev1[i] + ev1[i - 1]) * (sim_output.sample_time[i] - sim_output.sample_time[i - 1]) / 3600
            if ev1_energy < 0:
                print('ev1 finish at', sim_output.sample_time[i])


    plt.figure()
    plt.plot(sim_output.sample_time,total_power,'k', label='other loads')
    plt.plot(sim_output.sample_time,ev1,'r', label='EV_1')

    plt.plot(sim_output.sample_time, ev2,'c', label='EV_2')
    plt.plot(sim_output.sample_time, ev3,'b', label='EV_3')
    xlim = np.arange(sim_output.min_time, 60 * 60 * 24, 60 * 60 * 3)
    plt.ylabel("power (W)")
    plt.xlabel("hour")
    plt.xticks(xlim, [str(n).zfill(2) + ':00' for n in np.arange(int(sim_output.min_time / 3600), 24, 3)])

    plt.legend()
    plt.show()
    np.savetxt("1.csv", np.vstack((sim_output.sample_time, p2ce(sim_output.sample_time,ev1))).T, delimiter=' ', fmt="%d %f")
    np.savetxt("2.csv", np.vstack((sim_output.sample_time, p2ce(sim_output.sample_time,ev2))).T, delimiter=' ', fmt="%d %f")
    np.savetxt("3.csv", np.vstack((sim_output.sample_time, p2ce(sim_output.sample_time,ev3))).T, delimiter=' ', fmt="%d %f")
