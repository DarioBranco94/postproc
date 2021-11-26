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

import os
import csv
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
import json
from bin import utils as ut
from bin import ev as ev
from bin import node as nd
import time
import logging

logging.basicConfig(level=logging.WARNING)


class Checker:
    energyDict = {}
    astDict = {}
    num_of_timeseries = 0
    energyProducerDict = {}
    energyEVDict = {}
    totalEnergyConsumption = 0
    totalEnergyCharged = 0
    totalEnergyProduced = 0
    powerPeakListFiles = {}
    energyChargedWithIdAsKey = {}
    pvListResampled = {}
    consumerResampled = {}
    selfConsumedEnergy = 0
    totalProd = 0
    selfC = 0
    shareOfBatteryCapacity = 0
    peakLoadList = {}
    estlstList = {}
    batteryCapacity = {}
    listOfPeaks = {}
    root = '.'
    reachedLimits = {}
    cpNum = 0
    ast_lst_constraint = {}
    estlstList = {}
    energy_respected_to_capacity = {}
    energy_charged_respect_to_Connection = {}
    selfConsumedEnergyRespectToPVProduction = ''
    chargingPowerLowerThanMaxChPowConstraint = {}
    offeredFlexibilityIndex = 0
    actualFlexibilityIndex = 0
    V2GFlexibilityIndex = 0
    evList = {}
    selfsuff = 0
    peakToAverage = 0
    numcharging = 0
    KPI551 = 0
    KPI552 = 0
    KPI553 = 0
    KPI554 = 0
    KPI555 = 0
    KPI556 = 0
    meanPower = 0
    peakPower = 0
    peakTime = 0

    def doChecks(self, path, startTime, pathXML, pathVisualizer):
        """
        Args:
            path:
            startTime:
            pathXML:
            pathVisualizer:
        """
        try:
            logging.warning("removing " + path)
            os.remove(path + "/checks/outputParam.csv")
            os.remove(path + "/checks/kpi.csv")
        except:
            logging.warning("error removing " + path);

        allfiles = [os.path.join(dp, f) for dp, dn, filenames in os.walk(path) for f in filenames if f.endswith('.csv')]
        self.readConsumptionProduction(allfiles, self.energyDict)
        prod_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(path + "/PV") for f in filenames if
                      f.endswith('.csv')]
        self.readConsumptionProduction(prod_files, self.energyProducerDict)
        ev_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(path + "/EV") for f in filenames if
                    f.endswith('.csv')]
        self.readConsumptionProduction(ev_files, self.energyEVDict)
        for key, energy in self.energyDict.items():
            self.totalEnergyConsumption += float(energy)
        for key, energy in self.energyEVDict.items():
            self.totalEnergyCharged += float(energy)
        for key, energy in self.energyProducerDict.items():
            self.totalEnergyProduced += float(energy)
        self.totalEnergyConsumption = self.totalEnergyConsumption - self.totalEnergyProduced
        self.workWithOutputTXT(path)

        for key, value in self.powerPeakListFiles.items():
            if (value.split('/')[-2] != 'EV'):
                self.powerPeakListFiles[key] = ut.generatePowerTimeSeries(value, startTime)
            else:
                self.powerPeakListFiles[key] = ut.generatePowerTimeSeries(value, startTime)

        for key, energy in self.energyProducerDict.items():
            self.pvListResampled[key] = ut.generatePowerTimeSeries(key, startTime)
        for key, energy in self.energyDict.items():
            if (key not in self.pvListResampled):
                if (key.split('/')[-2] != 'EV'):
                    self.consumerResampled[key] = ut.generatePowerTimeSeries(key, startTime)
                else:
                    self.consumerResampled[key] = ut.generatePowerTimeSeries(key, startTime)
        self.calculateSelfConsumption()
        self.calculatePeakToAverage(startTime)
        self.readNeighborhoodXML(pathXML, startTime)

        self.readLoadXML(pathXML, startTime)
        # sumForPowerPeak(self.root, self.powerPeakListFiles)
        # findPeak(self.root, self.listOfPeaks)
        # self.checkPowerPeakConstraint()

        self.checkEnergyRespectToCapacityConstraint()
        self.checkEnergyRespectToConnectionTime()
        self.checkselfConsumedEnergyRespectToProduction()
        self.checkChargingPowerLowerThanMaxChPowConstraint(startTime)
        self.calculateSHareOfBatteryCapacityForV2G()
        self.calculateChargingFlexibility(startTime)
        self.UtilisationOfCps()
        self.calculateChargingAvailability(startTime)
        self.plotAll(startTime, path)
        self.writeOutput(path)

    def plotAll(self, startTime, path):
        endTime = startTime + 86400
        xnew = np.arange(startTime, endTime, 300)
        fig = plt.figure()
        xhh = []

        for element in xnew:
            xhh.append(time.strftime('%H:%M', time.localtime(element - 7200)))  # aggiusto la timezone
        for key, element in self.consumerResampled.items():
            plt.plot(xhh, element)
        for key, element in self.pvListResampled.items():
            plt.plot(xhh, element)
        plt.xticks(np.arange(0, 288, 35))
        fig.savefig(path + '/results.png', dpi=fig.dpi)

    def checkChargingPowerLowerThanMaxChPowConstraint(self, startTime):
        """
        Args:
            startTime:
        """
        for key, ev in self.evList.items():
            filename = ev.profile
            powerProfile = ut.generatePowerTimeSeries(filename, startTime)
            respected = 1
            for element in powerProfile:
                if (float(element) > float(ev.maxch)):
                    respected = 0
            if (respected == 0):
                self.chargingPowerLowerThanMaxChPowConstraint[key] = 'NotRespected'
            else:
                self.chargingPowerLowerThanMaxChPowConstraint[key] = 'Respected'

    def checkselfConsumedEnergyRespectToProduction(self):
        if (self.selfConsumedEnergy <= self.totalProd):
            self.selfConsumedEnergyRespectToPVProduction = 'Respected'
        else:
            self.selfConsumedEnergyRespectToPVProduction = 'Not Respected'

    def checkEnergyRespectToCapacityConstraint(self):
        for key, ev in self.evList.items():
            capacity = ev.capacity
            if (round(float(self.energyChargedWithIdAsKey[key]) + (float(ev.soc) * float(capacity)) / 100, 1) <= float(
                    capacity)):
                self.energy_respected_to_capacity[key] = 'Respected'
            else:
                self.energy_respected_to_capacity[key] = 'Not Respected'

    def checkEnergyRespectToConnectionTime(self):
        for key, ev in self.evList.items():
            capacity = ev.capacity
            if ((float(ev.adt) - (float(ev.aat))) * float(ev.maxch) > float(self.energyChargedWithIdAsKey[key])):
                self.energy_charged_respect_to_Connection[key] = 'Respected'
            else:
                self.energy_charged_respect_to_Connection[key] = 'Not Respected'

    def checkPowerPeakConstraint(self):
        for key, peak in self.peakLoadList.items():
            if (float(self.listOfPeaks[key]) * 1000 > float(peak)):
                self.reachedLimits[key] = 'reached'
            else:
                self.reachedLimits[key] = 'not reached'

    def checkASTConstraint(self):
        for key, ast in self.astDict.items():
            if (float(ast) < float(self.estlstList[key][1]) and float(ast) > float(self.estlstList[key][0])):
                self.ast_lst_constraint[key] = 'Respected'
            else:
                self.ast_lst_constraint[key] = 'not Respected'

    def readNeighborhoodXML(self, pathXML, startTime):
        """
        Args:
            pathXML:
            startTime:
        """
        tree = ET.parse(pathXML + '/neighborhood.xml')
        neighborhood = tree.getroot()
        for elem in neighborhood:  # CASE
            if (elem.tag != "fleet"):
                if (elem.tag == "chargingStation"):
                    self.cpNum += 1
                for subelement in elem:  # UTENTI
                    if (subelement.tag == "chargingStation"):
                        self.cpNum += 1
                    for subsubelement in subelement:  # charginPoint
                        if (subsubelement.tag == "chargingStation"):
                            self.cpNum += 1
            if (elem.tag == "fleet"):
                for ecar in elem:  # car
                    if (ecar.tag == "ecar"):
                        tempo = '[' + ecar.find("id").text + ']'
                        try:
                            self.evList[tempo].capacity = float(ecar.find("capacity").text)
                            self.evList[tempo].maxch = float(ecar.find("maxchpowac").text)
                            self.evList[tempo].minch = float(ecar.find("maxchpowac").text)  # DA CAMBIARE IN MIN
                            self.evList[tempo].maxDisPow = float(ecar.find("maxdispowac").text)

                        except:
                            None

    def readLoadXML(self, pathXML, startTime):
        """
        Args:
            pathXML:
            startTime:
        """
        tree = ET.parse(pathXML + '/loads.xml')
        neighborhood = tree.getroot()
        buildingID = "["
        self.root = nd.Node("root")
        for elem in neighborhood:
            buildingID = "["
            # print(elem)
            # elemNode = self.root.addChild("[" + elem.attrib['id'] + "]")
            if 'peakLoad' in elem.attrib:
                buildingID += elem.attrib['id'] + "]"
                self.peakLoadList[buildingID] = elem.attrib['peakLoad']
                buildingID += ":["
                for subelement in elem:
                    if 'peakLoad' in subelement.attrib:
                        tempo = buildingID + subelement.attrib['id'] + "]"
                        self.peakLoadList[tempo] = subelement.attrib['peakLoad']
                        # elemNode.addChild(tempo)
                    for subsubelement in subelement:
                        if 'peakLoad' in subsubelement.attrib:
                            tempo = buildingID + subsubelement.attrib['id'] + "]"
                            self.peakLoadList[tempo] = subsubelement.attrib['peakLoad']
                            # elemNode.addChild(tempo)
                        # elif (subsubelement.tag == "device"):
                        #    tempo = buildingID + subsubelement.find("id").text + ']'
                        #    if (subsubelement.find("est").text != "0" and subsubelement.find("lst").text != "0"):
                        #        self.estlstList[tempo] = [int(subsubelement.find("est").text) + int(startTime),
                        #                                  int(subsubelement.find("lst").text) + int(startTime)]
            for ecar in elem:
                if (ecar.tag == "ecar"):
                    self.numcharging += 1

                    tempo = '[' + ecar.find("id").text + ']'
                    self.evList[tempo].aat = float(ecar.find("aat").text)
                    self.evList[tempo].adt = float(ecar.find("adt").text)
                    self.evList[tempo].pat = float(ecar.find("pat").text)
                    self.evList[tempo].pdt = float(ecar.find("pdt").text)
                    self.evList[tempo].soc = float(ecar.find("soc").text)
                    self.evList[tempo].targetSoc = float(ecar.find("targetSoc").text)
                    self.evList[tempo].departureTimeMinusArrivalTimeMinusEnergyDemand = self.evList[tempo].maxch * (
                                (float(
                                    self.evList[tempo].adt) - self.evList[tempo].aat) / 3600 - (
                                            (self.evList[tempo].targetSoc - self.evList[tempo].soc) * self.evList[
                                        tempo].capacity / 100) / self.evList[tempo].maxch) * 0.5

    def calculatePeakToAverage(self, startTime):
        totalpowTS = []
        for i in range(288):
            tempCon = 0
            for key, power in self.consumerResampled.items():
                tempCon += power[i]
            totalpowTS.append(tempCon)
        energy = 0
        peakPow = totalpowTS[0]
        count = 0
        for element in totalpowTS:
            count += 1
            energy += element * (300 / 3600)
            if (element > peakPow):
                peakPow = element
                self.peakTime = count
        self.peakTime = self.peakTime * 300 + startTime
        self.peakTime = time.strftime('%H:%M', time.localtime(self.peakTime - 7200))
        meanPow = energy / 24
        self.meanPower = meanPow
        self.peakPower = peakPow
        if (meanPow != 0):
            self.peakToAverage = peakPow / meanPow
        else:
            self.peakToAverage = 9999999

    def calculateSelfConsumption(self):
        totalcon = 0
        for i in range(288):
            tempCon = 0
            tempProd = 0
            for key, power in self.consumerResampled.items():
                tempCon += power[i]
            for key, power in self.pvListResampled.items():
                tempProd += power[i]
            if (tempCon < tempProd):
                self.selfConsumedEnergy += tempCon
            else:
                self.selfConsumedEnergy += tempProd
            totalcon += tempCon
            self.totalProd += tempProd
        if (self.totalProd != 0):
            self.selfC = self.selfConsumedEnergy / self.totalProd
        else:
            self.selfC = 0
        if (totalcon != 0):
            self.selfsuff = self.selfConsumedEnergy / totalcon
        else:
            self.selfsuff = 0

    def readConsumptionProduction(self, allfiles, dictionary):
        """
        Args:
            allfiles:
            dictionary:
        """
        for file in allfiles:
            self.num_of_timeseries = self.num_of_timeseries + 1
            try:
                with open(file) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter=',')
                    for row in csv_reader:
                        if (row[1] != '0'):
                            dictionary[file.replace('\\', '/')] = row[1]
            except:
                None
        try:
            del dictionary["outputParam.csv"]
        except KeyError:
            pass

    def workWithOutputTXT(self, path):
        """
        Args:
            path:
        """
        file1 = open(path + '/output.txt', 'r')
        Lines = file1.readlines()
        for line in Lines:
            splittedMessage = line.split(" ")
            if (splittedMessage[0] == "<<<" or splittedMessage[0] == ">>>"):
                if (splittedMessage[1] == "ASSIGNED_START_TIME"):
                    self.astDict[splittedMessage[2]] = splittedMessage[4].rstrip()
                if (splittedMessage[1] == "LOAD"):
                    idList = splittedMessage[2].split(':')
                    idList.pop()
                    id = ''
                    first = 1
                    for value in idList:
                        if (first == 1):
                            id = value
                            first = 0
                        else:
                            id = id + ':' + value
                    csv_name = path + '/SH/' + splittedMessage[6].split('/')[-1]
                    self.powerPeakListFiles[id] = csv_name
                ev_dict = json.loads(line[4:])
                if (ev_dict['message']['subject'] == "EV"):
                    id = ev_dict['message']['id'].split(':')[1]
                    csv_name = path + '/EV/' + ev_dict['message']['id'].split(":")[0] + "_" + \
                               ev_dict['message']['id'].split(":")[1] + '.csv'
                    self.powerPeakListFiles[id] = csv_name
                    self.evList[id] = ev.EV()
                    self.energyChargedWithIdAsKey[id] = self.energyEVDict[csv_name]
                    self.evList[id].profile = csv_name

    def calculateChargingAvailability(self, startTime):
        chargedEnergy = 0
        chargeDemand = 0
        firstSample = 0
        lastSample = 0
        self.KPI551 = 0
        self.KPI552 = 0
        self.KPI553 = 0
        self.KPI554 = 0
        self.KPI555 = 0
        self.KPI556 = 0
        plugInDelay = 0
        notFinishedInTime = 0
        plugOutDelay = 0
        for key, ev in self.evList.items():
            with open(ev.profile) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                first = 0
                for row in csv_reader:
                    if (first == 0):
                        plug_in_time = float(row[0])
                        firstSample = float(row[1])
                        first = 1
                    plug_out_time = float(row[0])
                    lastSample = float(row[1])
            chargedEnergy += lastSample - firstSample
            if ((lastSample - firstSample) == 0):
                self.KPI553 += 1
            if ((ev.pat + startTime) > plug_in_time):
                plugInDelay += ev.pat + startTime - plug_in_time
            if (plug_out_time > ev.pdt + startTime):
                notFinishedInTime += 1
                plugOutDelay += plug_out_time - (ev.pdt + startTime)
            chargeDemand += ev.capacity * (ev.targetSoc - ev.soc) / 100
        self.KPI551 = chargedEnergy / chargeDemand
        if (chargedEnergy >= chargeDemand):
            self.KPI552 += 1
        self.KPI553 = self.KPI553 / self.numcharging
        self.KPI554 = plugInDelay / self.numcharging
        self.KPI555 = notFinishedInTime / self.numcharging
        self.KPI556 = plugOutDelay / self.numcharging / 3600

    def calculateSHareOfBatteryCapacityForV2G(self):
        for key, ev in self.evList.items():
            maxChPow = ev.maxch
            self.shareOfBatteryCapacity += 0.5 * maxChPow * ev.departureTimeMinusArrivalTimeMinusEnergyDemand
        if (self.shareOfBatteryCapacity <= 0.10):
            self.shareOfBatteryCapacity = 0

    def calculateChargingFlexibility(self, startTime):
        connectedTime = 0
        offeredFlexibility = []
        actualFlexibility = []
        V2GFlexibility = []
        for key, ev in self.evList.items():
            with open(ev.profile) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    plug_out_time = float(row[0])
            offeredFlexTime = (ev.adt - ev.aat) / 3600
            actualFlexTime = (plug_out_time - startTime - ev.aat) / 3600
            RequiredEnergy = ev.capacity * (ev.targetSoc - ev.soc) / 100

            if offeredFlexTime != 0:
                offeredFlexibility.append(1 - (RequiredEnergy / ev.maxch) / offeredFlexTime)
                V2GFlexibility.append(1 - (RequiredEnergy / ev.maxDisPow) / offeredFlexTime)
            else:
                offeredFlexibility.append(0)
                V2GFlexibility.append(0)
            if actualFlexTime != 0:
                actualFlexibility.append(1 - (RequiredEnergy / ev.maxch) / actualFlexTime)
            else:
                actualFlexibility.append(0)

        offeredFlexibilityIndex = 0
        actualFlexibilityIndex = 0
        V2GFlexibilityIndex = 0
        for i in range(0, len(offeredFlexibility)):
            offeredFlexibilityIndex += offeredFlexibility[i]
            actualFlexibilityIndex += actualFlexibility[i]
            V2GFlexibilityIndex += V2GFlexibility[i]
        if (len(offeredFlexibility) != 0):
            self.offeredFlexibilityIndex = offeredFlexibilityIndex / len(offeredFlexibility)
            self.actualFlexibilityIndex = actualFlexibilityIndex / len(offeredFlexibility)
            self.V2GFlexibilityIndex = V2GFlexibilityIndex / len(offeredFlexibility)
        else:
            self.offeredFlexibilityIndex = 0
            self.actualFlexibilityIndex = 0
            self.V2GFlexibilityIndex = 0
        if (self.offeredFlexibilityIndex <= 0):
            self.offeredFlexibilityIndex = 0
        if (self.actualFlexibilityIndex <= 0):
            self.actualFlexibilityIndex = 0
        if (self.V2GFlexibilityIndex <= 0):
            self.V2GFlexibilityIndex = 0

    def UtilisationOfCps(self):
        timespan = 86400
        connectedTime = 0
        chargingTime = 0
        chargedEnergy = 0
        for key, ev in self.evList.items():
            with open(ev.profile) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                first = 0
                for row in csv_reader:
                    if (first == 0):
                        plug_in_time = float(row[0])
                        first = 1
                    plug_out_time = float(row[0])
            connectedTime += ev.adt - ev.aat
            chargingTime += plug_out_time - plug_in_time
            chargedEnergy += ev.capacity * (ev.targetSoc - ev.soc)/100
        self.KPI531 = connectedTime / timespan
        self.KPI532 = chargingTime / connectedTime
        self.KPI533 = 3600*(chargedEnergy / connectedTime)

    def writeOutput(self, path):

        """
        Args:
            path:
        """
        try:
            os.mkdir(path + "/checks")
        except:
            pass
        '''
        with open(path+"/checks/outputParam.csv", "w") as csv_file:
                    param_writer = csv.writer(csv_file, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    param_writer.writerow(["Total_Energy_Consumption", str(totalEnergyConsumption)])
                    param_writer.writerow(["Total_Energy_Production", str(totalEnergyProduced)])
                    param_writer.writerow(["Assigned Start Time List", str(astDict)])
                    param_writer.writerow(["Number_of_Timeseries", str(num_of_timeseries)])
                    param_writer.writerow(["Energy_charged", str(totalEnergyCharged)])
                    param_writer.writerow(["Self_Consumption", str(selfC)])
        '''

        with open(path + "/checks/parameters.js", "w") as json_file:

            test_values = self.get_test_value()
            parameters = {"Total_Energy_Consumption": [str(self.totalEnergyConsumption), ""],
                          "Total_Energy_Production": [str(self.totalEnergyProduced), ""],
                          "Assigned Start Time List": [str(self.astDict), ""],
                          "Energy_charged": [str(self.totalEnergyCharged), ""],
                          "Number_of_Timeseries": [str(self.num_of_timeseries), ""],
                          "Self_Consumption": [str(self.selfC), ""]

                          }
            if test_values is not None:
                for test_value in test_values:
                    test_key = test_value[0].decode("utf8")
                    if test_key in parameters.keys():
                        parameters[test_key][1] = test_value[1].decode("utf8")

            json_file.write("data={rows:[")
            i = 1
            for key, value in parameters.items():
                json_file.write(
                    '{id:' + str(i) + ',data:[ "' + key + '","' + value[0] + '","' + value[1] + '","",""]},')
                i += 1
            json_file.write(']};')

        # temporary

        with open(path + "/checks/checks.js", "w") as json_file:
            json_file.write("coherence_checks={rows:[")
            json_file.write(
                '{id:1,data:[ "AstLstContraintRespected","' + str(self.ast_lst_constraint) + '","","",""]},')
            json_file.write('{id:2,data:[ "PowerPeaksReached","' + str(self.listOfPeaks) + '","","",""]},')
            json_file.write('{id:3,data:[ "PowerPeaksLimits","' + str(self.peakLoadList) + '","","",""]},')
            json_file.write('{id:4,data:[ "PowerPeaksLimitsReached","' + str(self.reachedLimits) + '","","",""]},')
            json_file.write('{id:5,data:[ "Energy_Charged_respect_to_capacity","' + str(
                self.energy_respected_to_capacity) + '","","",""]},')
            json_file.write('{id:6,data:[ "Energy_Charged_respect_to_Connection","' + str(
                self.energy_charged_respect_to_Connection) + '","","",""]},')
            json_file.write('{id:7,data:[ "Energy_AutoConsumed_Respect_To_Energy_Produced","' + str(
                self.selfConsumedEnergyRespectToPVProduction) + '","","",""]},')
            json_file.write('{id:8,data:[ "Charging_Power_Lower_than_Maximum","' + str(
                self.chargingPowerLowerThanMaxChPowConstraint) + '","","",""]}')
            json_file.write(']};')

        with open(path + "/checks/outputParam.csv", "w") as csv_file:
            param_writer = csv.writer(csv_file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            param_writer.writerow(["Total_Energy_Consumption", str(self.totalEnergyConsumption)])
            param_writer.writerow(["Total_Energy_Production", str(self.totalEnergyProduced)])
            param_writer.writerow(["Number_of_Timeseries", str(self.num_of_timeseries)])
            param_writer.writerow(["Energy_charged", str(self.totalEnergyCharged)])
            param_writer.writerow(["Self_Consumption", str(self.selfC)])
            param_writer.writerow(["Assigned Start Time List", str(self.astDict)])
            param_writer.writerow(["ast_lst_List", str(self.estlstList)])
            param_writer.writerow(["AstLstContraintRespected", str(self.ast_lst_constraint)])
            param_writer.writerow(["PowerPeaksReached", str(self.listOfPeaks)])
            param_writer.writerow(["PowerPeaksLimits", str(self.peakLoadList)])
            param_writer.writerow(["PowerPeaksLimitsReached", str(self.reachedLimits)])
            param_writer.writerow(["Energy_Charged_respect_to_capacity", str(self.energy_respected_to_capacity)])
            param_writer.writerow(
                ["Energy_Charged_respect_to_Connection", str(self.energy_charged_respect_to_Connection)])
            param_writer.writerow(
                ["Energy_AutoConsumed_Respect_To_Energy_Produced", str(self.selfConsumedEnergyRespectToPVProduction)])
            param_writer.writerow(
                ["Charging_Power_Lower_than_Maximum", str(self.chargingPowerLowerThanMaxChPowConstraint)])

        with open(path + "/checks/kpi.csv", "w") as csv_file:
            param_writer = csv.writer(csv_file, delimiter=' ', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            param_writer.writerow(["NumberOfEV_GC5.1", str(len(self.evList))])
            param_writer.writerow(["NumberOfCPGC5.2", str(self.cpNum)])
            param_writer.writerow(["Self_Consumption_GC5.14.1", str(self.selfC * 100)])
            param_writer.writerow(["Self_Consumption_GC5.14.2", str(self.selfsuff * 100)])
            param_writer.writerow(["UtilisationOfCPsGC5.3.1", str(self.KPI531)])
            param_writer.writerow(["UtilisationOfCPsGC5.3.2", str(self.KPI532)])
            param_writer.writerow(["UtilisationOfCPsGC5.3.3", str(self.KPI533)])
            param_writer.writerow(["ChargingFlexibility5.13.1", str(self.offeredFlexibilityIndex)])
            param_writer.writerow(["ChargingFlexibility5.13.2", str(self.actualFlexibilityIndex)])
            param_writer.writerow(["ChargingFlexibility5.13.3", str(self.V2GFlexibilityIndex)])
            param_writer.writerow(["shareOfBatteryCapacity5.4", str(self.shareOfBatteryCapacity)])
            param_writer.writerow(["PeakToAverage5.10.1", str(self.peakToAverage)])
            param_writer.writerow(["PeakToAverage (AVERAGE POWER) 5.10.2", str(self.meanPower)])
            param_writer.writerow(["PeakToAverage (PEAK POWER) 5.10.3", str(self.peakPower)])
            param_writer.writerow(["PeakToAverage (PEAK Time) 5.10.4", str(self.peakTime)])
            param_writer.writerow(["NumOfCharging5.3.4", str(self.numcharging)])
            param_writer.writerow(["ChargingAvailability5.5.1", str(self.KPI551)])
            param_writer.writerow(["ChargingAvailability5.5.2", str(self.KPI552)])
            param_writer.writerow(["ChargingAvailability5.5.3", str(self.KPI553)])
            param_writer.writerow(["ChargingAvailability5.5.4", str(self.KPI554)])
            param_writer.writerow(["ChargingAvailability5.5.5", str(self.KPI555)])
            param_writer.writerow(["ChargingAvailability5.5.6", str(self.KPI556)])
        with open(path + "/checks/kpis.js", "w") as json_file:
            json_file.write("kpis_values={rows:[")
            json_file.write('{id:1,data:[ "GC5.1","Number Of EVs",' + str(len(self.evList)) + ']},')
            json_file.write('{id:2,data:[ "GC5.2","Number Of CP",' + str(self.cpNum) + ']},')
            json_file.write('{id:3,data:[ "GC5.14.1","Self Consumption",' + str(self.selfC * 100) + ']},')
            json_file.write('{id:4,data:[ "GC5.14.2","Self Consumption",' + str(self.selfsuff * 100) + ']},')

            json_file.write('{id:5,data:[ "GC5.3.1","Utilisation Of CPs",' + str(self.KPI531) + ']},')
            json_file.write('{id:6,data:[ "GC5.3.2","Utilisation Of CPs",' + str(self.KPI532) + ']},')
            json_file.write('{id:7,data:[ "GC5.3.3","Utilisation Of CPs",' + str(self.KPI533) + ']},')
            json_file.write(
                '{id:8,data:[ "GC5.13.1","Charging Flexibility",' + str(self.offeredFlexibilityIndex) + ']},')
            json_file.write(
                '{id:9,data:[ "GC5.13.2","Charging Flexibility",' + str(self.actualFlexibilityIndex) + ']},')
            json_file.write('{id:9,data:[ "GC5.13.3","Charging Flexibility",' + str(self.V2GFlexibilityIndex) + ']},')
            json_file.write(
                '{id:10,data:[ "GC5.4","Share Of Battery Capacity",' + str(self.shareOfBatteryCapacity) + ']},')
            json_file.write(
                '{id:11,data:[ "GC5.10","Peak To Average",' + str(self.peakToAverage) + ']},')
            json_file.write(
                '{id:12,data:[ "GC5.3.4","Num of Charging",' + str(self.numcharging) + ']},')
            json_file.write('{id:13,data:[ "GC5.5.1","Charging Availability",' + str(self.KPI551) + ']},')
            json_file.write('{id:14,data:[ "GC5.5.2","Charging Availability",' + str(self.KPI552) + ']},')
            json_file.write('{id:15,data:[ "GC5.5.3","Charging Availability",' + str(self.KPI553) + ']},')
            json_file.write('{id:16,data:[ "GC5.5.4","Charging Availability",' + str(self.KPI554) + ']},')
            json_file.write('{id:17,data:[ "GC5.5.5","Charging Availability",' + str(self.KPI555) + ']},')
            json_file.write('{id:18,data:[ "GC5.5.6","Charging Availability",' + str(self.KPI556) + ']},')
            json_file.write('{id:19,data:[ "GC5.10.2","Average Power",' + str(self.meanPower) + ']},')
            json_file.write('{id:20,data:[ "GC5.10.3","Peak Power",' + str(self.peakPower) + ']},')
            json_file.write('{id:21,data:[ "GC5.10.4","Peak Time",' + str(self.peakTime.split(':')[0]+'.' + self.peakTime.split(':')[1]) + ']}')
            json_file.write(']};')
    def get_test_value(self):
        cwd = os.getcwd()
        cwd_parts = cwd.split("/")
        sim_date = cwd_parts[-1]
        parts = sim_date.split("_")

        sim_date = ' 12_12_15'
        test_file = "../../tests/" + sim_date + "/outputParam.csv"
        test_values = None
        if os.path.isfile(test_file):
            test_values = np.genfromtxt(test_file, delimiter=";", dtype="|S")
        return test_values







