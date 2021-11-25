
from os import path
import os
import csv
import glob
import xml.etree.ElementTree as ET
import shutil
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
from bin import visualization
import sqlite3
from bin import config
import datetime
import json



def findPeak(node, listOfPeaks):
    """
    Args:
        node:
        listOfPeaks:
    """

    max = 0
    for i in range(0, 287):
        if (float(node.data[i][0]) > max):
            max = node.data[i][0]
    if (node.name != 'root'):
        listOfPeaks[node.name] = max

    for nodechild in node.children:
        findPeak(nodechild, listOfPeaks)


def printChilds(node):
    """
    Args:
        node:
    """

    for i in range(0, 287):
        print(node.data[i])
    for childNode in node.children:
        printChilds(childNode)


def sumForPowerPeak(node, dictConsumer):
    """
    Args:
        node:
        dictConsumer:
    """
    for nodechild in node.children:

        # print(nodechild.name)
        powerList = sumForPowerPeak(nodechild, dictConsumer)
        for i in range(0, 287):
            node.data[i] += powerList[i]

    for key, consumer in dictConsumer.items():
        if (key == node.name):
            # print(key)
            for i in range(0, 287):
                node.data[i] += consumer[i]
    return node.data


def generateEnergyTimeSeries(file, startTime):
    """
    Args:
        file:
        startTime:
    """
    endTime = startTime + 86400
    with open(file, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=",")
        count = 0
        x = []  # lista dei tempi della timeseries
        y = []  # lista dei valori della timeseries
        lastSample = 0
        lastValue = 0  # Questo mi serve per tenermi in memoria il valore di energia precedente alla riga che sto leggendo, cosi posso farmi il delta per la trasformazione in potenza
        for row in csv_reader:  # per tutte le righe
            if (count != 0):  # salto la prima riga della ts
                if (float(row[1]) != 0):
                    x.append(float(row[0]))  # aggiunto il tempo alla lista dei tempi
                    y.append((float(row[1]) - lastValue))
                else:
                    x.append(float(row[0]))  # aggiunto il tempo alla lista dei tempi
                    y.append((float(row[1])))
            else:
                if (startTime < float(row[
                                          0])):  # faccio in modo che se il primo tempo della timeseries Ã© piÃº grande del minimo del periodo di interesse ci piazzo uno zero, cosi dopo non ho problemi quando vado a ricampionare
                    x.append(startTime)
                    y.append(0)
                else:
                    x.append(float(row[0]))  # aggiunto il tempo alla lista dei tempi
                    y.append(float(row[1]))  # aggiungo alla lista dei valori la potenza
            lastSample = float(row[0])
            lastValue = float(row[1])  # aggiorno l'energia precedente
            count += 1  # aggiorno il count quando ho finito la riga
    if (
            endTime > lastSample):  # stesso discorso di prima, se l'ultimo tempo della timeseries Ã© piÃº piccolo del massimo tempo di interesse metto uno zero per non aver problemi dopo
        y.append(0)
        x.append(endTime)
    f = interpolate.interp1d(x, y)  # faccio l'interpolazione lineare
    xnew = np.arange(startTime, endTime, 300)  # mi creo il vettore dei tempi con un sample ogni 5 minuti (300 secondi)

    ynew = f(xnew)

    return ynew


def generatePowerTimeSeries(file, startTime):
    """
    Args:
        file:
        startTime:
    """

    endTime = startTime + 86400
    with open(file, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=",")
        count = 0

        x = []  # lista dei tempi della timeseries
        y = []  # lista dei valori della timeseries
        lastSample = 0  # Questo mi serve per tenermi in memoria il tempo precedente alla riga che sto leggendo, cosi posso farmi il delta per la trasformazione in potenza
        lastValue = 0  # Questo mi serve per tenermi in memoria il valore di energia precedente alla riga che sto leggendo, cosi posso farmi il delta per la trasformazione in potenza
        for row in csv_reader:  # per tutte le righe
            try:
                row[0] = datetime.datetime.strptime( row[0], "%Y-%m-%d %H:%M:%S").timestamp()
            except:
                None
            if (count != 0):  # salto la prima riga del
                # la ts perchÃ© devo convertire in potenza
                if(float(row[0]) != lastSample):
                    x.append(lastSample+1)  # aggiunto il tempo alla lista dei tempi
                    y.append(3600*(float(row[1]) - lastValue) / (float(row[0]) - lastSample))
                    x.append(float(row[0]))  # aggiunto il tempo alla lista dei tempi
                    y.append(3600 * (float(row[1]) - lastValue) / (float(row[0]) - lastSample))

            else:
                if (startTime < float(row[0])):  # faccio in modo che se il primo tempo della timeseries Ã© piÃº grande del minimo del periodo di interesse ci piazzo uno zero, cosi dopo non ho problemi quando vado a ricampionare
                    x.append(startTime)
                    y.append(0)  # aggiungo alla lista dei valori la potenza
            lastSample = float(row[0])  # aggiorno il tempo precedente
            lastValue = float(row[1])  # aggiorno l'energia precedente
            count += 1  # aggiorno il count quando ho finito la riga
    y.append(0)
    x.append(lastSample+1)
    #if endTime > lastSample:  # stesso discorso di prima, se l'ultimo tempo della timeseries Ã© piÃº piccolo del massimo tempo di interesse metto uno zero per non aver problemi dopo
    y.append(0.1)
    x.append(endTime)
    f = interpolate.interp1d(x, y, kind='previous')  # faccio l'interpolazione lineare
    xnew = np.arange(startTime, endTime, 300)  # mi creo il vettore dei tempi con un sample ogni 5 minuti (300 secondi)
    ynew = f(xnew)  # genero la nuova serie di potenze ricampionatew
    #plt.plot(xnew,ynew)
    #plt.show()
    return ynew


def html_images(folder):
    """
    Args:
        folder:
    """
    html_file = open(folder + "/index.html", "w")
    html_file.write('<html><body>')
    images = glob.glob("./*.png")
    for image in images:
        shutil.copy2(image, folder)
        os.remove(os.path.basename(image))
    images = glob.glob(folder + "/*.png")
    for image in images:
        html_file.write('<img src="' + os.path.basename(image) + '"/>')
    html_file.write('</body></html>')
    html_file.close()

