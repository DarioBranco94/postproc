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

"""
Config
=======================================
This class is used to load the configuration file.
"""


import yaml
import os



######################################################################################
# This class manages scenario's configuration and static data written in config.yml. #
######################################################################################
class Configuration:
    """
    Configuration loader class.
    """
    listDevice = []  # IN THIS LIST WILL BE STORED ALL THE LOADS
    listPanels = []  # IN THIS LIST WILL BE STORED ALL THE PRODUCERS
    listEvent = []
    pathneighbor = 0
    pathload2 = 0
    print(os.listdir('../'))
    config_file = './config.yml'
    parameters = None
    mydir = None
    messageToWait = None

    @classmethod
    def set_config_file(cls, config_file='config.yml'):
        """
        Set the name of the configuration file to be read.
        Args:
            config_file: The path of configuration file.
        """
        cls.config_file = config_file

    @classmethod
    def load(cls):
        """
        This method loads the configuration file.
        """
        if cls.config_file is None:
            raise Exception("Config File is not set")
        with open(cls.config_file, 'r') as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.Loader)
            cls.parameters = cfg['config']
            cls.parameters['current_sim_dir'] = cls.parameters['simulation_dir'] + '/' + cls.parameters['simulation']
            cls.parameters['webdir'] = cfg['config']['webdir']
            if(cls.parameters['forcedwait'] == False):
                cls.messageToWait = ['load']
            else:
                cls.messageToWait = ['load', 'background', 'EV']

            mydir = 0
            # codice aggiunto perchÃ¨ non funzionava la directory web nel dispatcher
            date3 = cls.parameters['date'] .split()
            newdir2 = date3[0].replace('/','_')
            sim_temp2 = newdir2.split("_")
            lock1 = False
            if len(sim_temp2[0]) == 1:
                sim_temp2[0] =  "0" + sim_temp2[0]
                lock1 = True
            if len(sim_temp2[1]) == 1:
                sim_temp2[1] = "0" + sim_temp2[1]
                lock1 = True
            if lock1:
                newdir2 = sim_temp2[0] + "_" + sim_temp2[1] + "_" + sim_temp2[2]
            dirCount1 = 1
            while os.path.exists(cls.parameters['current_sim_dir'] + "/Results/" + newdir2 + "_" + str(dirCount1)):
                dirCount1 += 1
            cls.parameters['runtime_dir'] = cls.parameters['current_sim_dir'] + "/Results/" + newdir2+"_"+str(dirCount1)
            temp =  cls.parameters['runtime_dir'].split("/")
            cls.parameters['user_dir'] = temp[-1]
            cls.dirCount1 = dirCount1
            # fine codice aggiunto









