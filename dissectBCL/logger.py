import logging

log = logging.getLogger()

# Definitions.
def setLog(logFile):
    logging.basicConfig(
        filename=logFile,
        level="DEBUG",
        format="%(levelname)s    %(asctime)s    %(message)s",
        filemode='w'
    )
    log = logging.getLogger()
    log.info("Log Initiated.")