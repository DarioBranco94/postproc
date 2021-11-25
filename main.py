from bin import postproc as pp, utils as ut
import yaml

if __name__ == "__main__":
    checker = pp.Checker()
    config_file = open("config.yml")
    config_parameters = yaml.load(config_file, Loader=yaml.FullLoader)
    #visualization.callExternal("./08_16_21_11/")
    checker.doChecks("./Simulations/"+config_parameters['simName']+"/output", 1629072000,"./Simulations/"+config_parameters['simName']+"/xml","./08_16_21_11")
    ut.html_images("./")
