import pandas as pd


def parseUnclassifiedPct(reportPath):
    """
    Returns the % of reads reported as 'unclassified' in a kraken2
    --report TSV, or None if the report is empty/unreadable, or has
    no unclassified row.
    """
    try:
        reportDF = pd.read_csv(reportPath, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return None
    unclassified = reportDF[reportDF[5].str.strip() == "unclassified"]
    if unclassified.empty:
        return None
    return float(unclassified.iloc[0][0])


def pickThreshold(libraryType, config):
    """
    Returns the unclassified-% threshold (a float) that applies to
    libraryType: config['screening']['relaxed_threshold'] if libraryType
    is listed in config['screening']['relaxed_library_types'] (matched
    case-insensitively), else config['screening']['unclassified_threshold'].
    """
    relaxedTypes = {
        t.strip().lower()
        for t in config["screening"]["relaxed_library_types"].split(",")
        if t.strip()
    }
    if isinstance(libraryType, str) and libraryType.strip().lower() in relaxedTypes:
        return config["screening"].getfloat("relaxed_threshold")
    return config["screening"].getfloat("unclassified_threshold")


def needsEscalation(reportPath, libraryType, config):
    """
    Returns True if reportPath's unclassified % exceeds the threshold
    that applies to libraryType (see pickThreshold).
    """
    pct = parseUnclassifiedPct(reportPath)
    if pct is None:
        return False
    return pct > pickThreshold(libraryType, config)
