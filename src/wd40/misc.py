import configparser


def getConf(configfile):
    config = configparser.ConfigParser()
    config.read(configfile)
    return(dict(config))
