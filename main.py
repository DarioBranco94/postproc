from bin import postproc as pp, utils as ut

if __name__ == "__main__":
    checker = pp.Checker()

    mydir = "AS_addPPs_04-10-21_16-57-06"

    #visualization.callExternal("./08_16_21_11/")

    checker.doChecks("./Simulations/08_16_21_11_booking/output", 1629072000,"./Simulations/08_16_21_11_booking/xml","./08_16_21_11")
    #shutil.copy('../../../../../../dockers/gcsim/gcsimulator/templates/checks.html', 'checks.html')
    ut.html_images("./")
